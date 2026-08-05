// Shared with the main thread's WebGPU Patterson map path - see symmetry_utils.js
// for normalizeSGSymbol(), findSpaceGroupSetting(), and expandReflections().
importScripts('symmetry_utils.js');

// Set by calculatePattersonMap() so the message handler can report which
// operator set the expansion actually used, and any warnings raised.
let lastExpansionSymmetry = null;

/**
 * Radius, in Angstroms, around u=v=w=0 that is treated as the Patterson origin
 * peak. The origin is a genuine feature (its height is the sum of all
 * intensities) but it is one to two orders of magnitude above every
 * interatomic vector, so including it in min/max statistics collapses every
 * real peak into the bottom few percent of the scale. It is masked out of
 * statistics and of peak searching, but the map values themselves are left
 * untouched.
 */
const ORIGIN_MASK_ANGSTROM = 1.1;

/**
 * Builds a boolean mask, one entry per voxel, true where the voxel lies within
 * ORIGIN_MASK_ANGSTROM of the origin (periodic images included).
 */
function buildOriginMask(res, cell) {
    const orth = sharkoOrthMatrix(cell);
    const mask = new Uint8Array(res * res * res);
    for (let iw = 0; iw < res; iw++) {
        for (let iv = 0; iv < res; iv++) {
            for (let iu = 0; iu < res; iu++) {
                const d = sharkoFracToCartLength(iu / res, iv / res, iw / res, orth);
                if (d <= ORIGIN_MASK_ANGSTROM) mask[iw * res * res + iv * res + iu] = 1;
            }
        }
    }
    return mask;
}

/** min/max of a map ignoring masked voxels and non-finite values. */
function maskedExtrema(map, mask) {
    let maxVal = -Infinity, minVal = Infinity, n = 0;
    for (let i = 0; i < map.length; i++) {
        if (mask && mask[i]) continue;
        const v = map[i];
        if (!isFinite(v)) continue;
        if (v > maxVal) maxVal = v;
        if (v < minVal) minVal = v;
        n++;
    }
    return { minVal, maxVal, count: n };
}

/**
 * Solves a coordinate string like "u/2" based on a peak's (u,v,w).
 */
// Compiled solvers are cached by expression: the same handful of strings are
// evaluated once per (peak x Harker section), so with 50 peaks and a dozen
// sections the old code was calling the Function constructor - a full parse
// and compile - some 600 times per run to build the same few functions.
const _solverCache = new Map();
function getSolver(sanitizedExpression) {
    let fn = _solverCache.get(sanitizedExpression);
    if (!fn) {
        fn = new Function('u', 'v', 'w', `return ${sanitizedExpression}`);
        _solverCache.set(sanitizedExpression, fn);
    }
    return fn;
}

function solveCoordinate(solverString, peak) {
    if (solverString === '?') return '?';
    try {
        let sanitizedExpression = solverString.replace(/[^uvw\d\+\-\*\/\%\.\(\)\s]/g, '');
        // Convert double negatives to a plus sign
        sanitizedExpression = sanitizedExpression.replace(/--/g, '+'); 
        
        const solverFunc = getSolver(sanitizedExpression);
        const result = solverFunc(peak.u, peak.v, peak.w);
        if (typeof result !== 'number' || isNaN(result) || !isFinite(result)) { throw new Error(`Solver returned non-finite: ${result}`); }
        return (((result % 1) + 1) % 1).toFixed(3);
    } catch (e) { console.error(`Error solving: "${solverString}" for peak (${peak.u.toFixed(3)}, ${peak.v.toFixed(3)}, ${peak.w.toFixed(3)}):`, e); return 'err'; }
}

/**
 * Calculates the 3D Patterson map by FFT.
 *
 * This used to be a direct Fourier summation: a triple loop over voxels with
 * an inner loop over every expanded reflection, so O(res^3 * numReflections).
 * At the default 50^3 grid with a typical 30000-reflection expanded list that
 * is nearly nine billion cosine evaluations, which is why it took minutes and
 * why a WebGPU version of the same algorithm was bolted on beside it.
 *
 * The summation is a discrete Fourier transform, so it can be done as one:
 * O(res^3 log res), around 340x faster measured on a 64^3 grid - fast enough
 * on a single worker thread that the GPU path is no longer worth its
 * complexity and has been removed.
 *
 * Returns { map, res, dMin, sigma } because the FFT chooses its own grid size
 * (a power of two large enough to avoid aliasing the highest-order
 * reflection), which is not necessarily the resolution that was requested.
 */
function calculatePattersonMap(crystalData, spaceGroups, mapResolution) {
    try {
        const { cell, reflections, spaceGroup } = crystalData;
        if (!reflections || reflections.length === 0) { throw new Error("No reflection data."); }
        if (!cell || !cell.a || !cell.b || !cell.c || isNaN(cell.a) || isNaN(cell.b) || isNaN(cell.c)) { throw new Error("Invalid cell data."); }

        // Expand unique reflections to the full sphere. Throws (rather than
        // silently expanding with the identity) if the symmetry database has
        // no usable operators for this setting.
        const expansion = expandReflections(reflections, spaceGroup.number, spaceGroup.name, spaceGroups,
                                            { perReflection: !!crystalData.perReflection });
        const fullReflections = expansion.reflections;
        lastExpansionSymmetry = expansion.symmetry;

        // General triclinic volume, shared with the rest of the program. The
        // old a*b*c is only correct for orthogonal cells; it is a pure scale
        // factor on the map, but the swarm's anti-bump penalty is scaled by
        // the map's own value range, so getting it wrong desynchronises
        // fitness from penalty.
        const result = sharkoPattersonFFT(fullReflections, cell, mapResolution);

        // Grid-size changes and dropped reflections are surfaced through the
        // same channel as the symmetry warnings so they cannot pass silently.
        if (result.warnings.length && lastExpansionSymmetry) {
            lastExpansionSymmetry.warnings = (lastExpansionSymmetry.warnings || []).concat(result.warnings);
        }
        return result;
    } catch (error) {
        console.error("[Worker] Map calc error:", error);
        throw error; // Re-throw to be caught by the main handler
    }
}

/**
 * Finds peaks in the calculated 3D map.
 */
function findPeaks(pattersonMap3D, mapResolution, cell) {
    const res = mapResolution, map = pattersonMap3D;
    if (!map) { return []; }

    // The Patterson function is periodic, so u=0, v=0 and w=0 are not edges -
    // they are the middle of the function, and they carry the Harker lines and
    // planes at 0 that this program exists to search. Scanning 1..res-2 made
    // every peak on those three faces unfindable by construction. Neighbour
    // indices wrap instead.
    const mask = cell ? buildOriginMask(res, cell) : null;
    const { minVal, maxVal } = maskedExtrema(map, mask);
    if (!isFinite(maxVal) || !isFinite(minVal) || maxVal === minVal) {
        console.warn(`[Worker] Map flat/invalid outside the origin. Skipping peaks.`);
        return [];
    }

    // Threshold and normalisation are now set by the interatomic vectors, not
    // by the origin peak (which is the sum of all intensities and would push
    // the 15% cut above every real feature).
    const threshold = minVal + (maxVal - minVal) * 0.15;
    const peaks = [];

    for (let iw = 0; iw < res; iw++) {
        for (let iv = 0; iv < res; iv++) {
            for (let iu = 0; iu < res; iu++) {
                const idx = iw * res * res + iv * res + iu;
                if (mask && mask[idx]) continue;          // inside the origin peak
                const val = map[idx];
                if (!isFinite(val) || val < threshold) continue;

                // A strict `nv > val` test makes every voxel of a flat plateau
                // a maximum, so a broad peak sitting on a few equal-valued
                // voxels was reported as several separate peaks that then
                // crowded out real ones in the 50-peak cut. Ties are broken by
                // index so exactly one voxel of any plateau survives.
                let isMax = true;
                for (let dw = -1; dw <= 1 && isMax; dw++) {
                    const jw = (iw + dw + res) % res;
                    for (let dv = -1; dv <= 1 && isMax; dv++) {
                        const jv = (iv + dv + res) % res;
                        for (let du = -1; du <= 1 && isMax; du++) {
                            if (du === 0 && dv === 0 && dw === 0) continue;
                            const ju = (iu + du + res) % res;
                            const nIdx = jw * res * res + jv * res + ju;
                            const nv = map[nIdx];
                            if (!isFinite(nv)) continue;
                            if (nv > val || (nv === val && nIdx < idx)) { isMax = false; break; }
                        }
                    }
                }
                if (isMax) {
                    peaks.push({
                        u: iu / res, v: iv / res, w: iw / res,
                        height: (val - minVal) / (maxVal - minVal),
                        value: val
                    });
                }
            }
        }
    }

    peaks.sort((a, b) => b.height - a.height);
    const foundPeaks = peaks.slice(0, 50);
    console.log(`[Worker] Found ${peaks.length} peaks (origin masked to ${ORIGIN_MASK_ANGSTROM} A). Kept ${foundPeaks.length}.`);
    return foundPeaks;
}

/**
 * Checks found peaks against Harker sections.
 */
function analyzeHarkerPeaks(foundPeaks, crystalData, spaceGroups, mapResolution) {
    let harkerAnalysisResults = [];
    if (!crystalData?.spaceGroup || foundPeaks.length === 0) { console.log("[Worker] Skipping Harker."); return []; }
    const sgNumber = crystalData.spaceGroup.number; 
    
    // Prioritize embedded sections from the Pawley file
    let sections = [];
    if (crystalData.harkerSections && crystalData.harkerSections.length > 0) {
        sections = crystalData.harkerSections;
        console.log(`[Worker] Using ${sections.length} embedded Harker sections.`);
    } else {
        // harker_sections lives inside the matched settings[] entry, not at
        // the top level of spaceGroups[sgNumber].
        const setting = findSpaceGroupSetting(spaceGroups, sgNumber, crystalData.spaceGroup.name);
        if (setting && setting.harker_sections && setting.harker_sections.length > 0) {
            sections = setting.harker_sections;
            console.log(`[Worker] Using JSON Harker sections for SG ${sgNumber} (${setting.symbol || 'setting'}).`);
        } else {
            console.warn(`[Worker] No Harker data available.`);
            return [];
        }
    }
    
    const gridSpacing = 1.0 / mapResolution; const tol = 1.5 * gridSpacing;
    console.log(`[Worker] Analyzing SG: ${sgNumber}. Tol: ${tol.toFixed(3)}`);
    
    sections.forEach((section, si) => {
        if (!section.coordinate || !['u', 'v', 'w'].includes(section.coordinate) || typeof section.value !== 'number' || !section.solver) { console.warn(`[Worker] Skip invalid section ${si + 1}`); return; }
        foundPeaks.forEach((peak, pi) => {
            const pc = peak[section.coordinate]; const diff = Math.abs(pc - section.value); const pDiff = Math.min(diff, 1.0 - diff);
            if (pDiff < tol) {
                const site = { source: `${section.type?.charAt(0).toUpperCase() + section.type?.slice(1) || 'Unk'} (${section.coordinate}=${section.value.toFixed(3)})`, peakCoords: `(${peak.u.toFixed(3)}, ${peak.v.toFixed(3)}, ${peak.w.toFixed(3)})`, x: solveCoordinate(section.solver.x, peak), y: solveCoordinate(section.solver.y, peak), z: solveCoordinate(section.solver.z, peak) };
                if (site.x === 'err' || site.y === 'err' || site.z === 'err') { console.error(`[Worker]   Solver error. Peak ${pi}, Sec ${si + 1}. Discarded.`); }
                else { harkerAnalysisResults.push(site); }
            }
        });
    });
    console.log(`[Worker] Harker found ${harkerAnalysisResults.length} partial site(s).`);
    return harkerAnalysisResults;
}


// --- Site Combination Helpers (for worker) ---
function averagePeriodic(v1, v2) { const diff = v1 - v2; if (Math.abs(diff) > 0.5) { if (v1 < v2) v1 += 1.0; else v2 += 1.0; } return ((( (v1 + v2) / 2.0 ) % 1) + 1) % 1; }
function adjustPeriodic(value, ref) { if (value - ref > 0.5) return value - 1.0; if (ref - value > 0.5) return value + 1.0; return value; }

/**
 * Combines partial Harker sites into full 3D atom sites.
 */
function combineSites(harkerAnalysisResults, harkerTolerance) {
    console.log("[Worker] --- Starting Site Combination ---");
    let consolidatedSites = [];
    const results = harkerAnalysisResults.filter(site => site.x !== 'err' && site.y !== 'err' && site.z !== 'err');
    const tol = harkerTolerance;
    console.log(`[Worker] Attempting to combine ${results.length} valid partial sites. Tolerance: ${tol.toFixed(3)}`);

    if (results.length < 2) {
        console.log("[Worker] --- Finished Site Combination (Not enough sites) ---");
        return [];
    }

    const areClose = (c1, c2) => { if (c1 === '?' || c2 === '?') return false; const v1 = parseFloat(c1); const v2 = parseFloat(c2); if (isNaN(v1) || isNaN(v2)) return false; const diff = Math.abs(v1 - v2); return Math.min(diff, 1 - diff) < tol; };
    const isNum = (c) => c !== '?' && !isNaN(parseFloat(c));

    const potentialSites = [];
    for (let i = 0; i < results.length; i++) {
        for (let j = i + 1; j < results.length; j++) {
            const r1 = results[i], r2 = results[j]; let combinedSite = null;
            try {
                if (areClose(r1.z, r2.z) && isNum(r1.x) && isNum(r2.y)) { const avgZ = averagePeriodic(parseFloat(r1.z), parseFloat(r2.z)); combinedSite = { x: parseFloat(r1.x), y: parseFloat(r2.y), z: avgZ }; }
                else if (areClose(r1.z, r2.z) && isNum(r2.x) && isNum(r1.y)) { const avgZ = averagePeriodic(parseFloat(r1.z), parseFloat(r2.z)); combinedSite = { x: parseFloat(r2.x), y: parseFloat(r1.y), z: avgZ }; }
                else if (areClose(r1.y, r2.y) && isNum(r1.x) && isNum(r2.z)) { const avgY = averagePeriodic(parseFloat(r1.y), parseFloat(r2.y)); combinedSite = { x: parseFloat(r1.x), y: avgY, z: parseFloat(r2.z) }; }
                else if (areClose(r1.y, r2.y) && isNum(r2.x) && isNum(r1.z)) { const avgY = averagePeriodic(parseFloat(r1.y), parseFloat(r2.y)); combinedSite = { x: parseFloat(r2.x), y: avgY, z: parseFloat(r1.z) }; }
                else if (areClose(r1.x, r2.x) && isNum(r1.y) && isNum(r2.z)) { const avgX = averagePeriodic(parseFloat(r1.x), parseFloat(r2.x)); combinedSite = { x: avgX, y: parseFloat(r1.y), z: parseFloat(r2.z) }; }
                else if (areClose(r1.x, r2.x) && isNum(r2.y) && isNum(r1.z)) { const avgX = averagePeriodic(parseFloat(r1.x), parseFloat(r2.x)); combinedSite = { x: avgX, y: parseFloat(r2.y), z: parseFloat(r1.z) }; }
                if (combinedSite) { const norm = val => (((val % 1) + 1) % 1); combinedSite.x = norm(combinedSite.x); combinedSite.y = norm(combinedSite.y); combinedSite.z = norm(combinedSite.z); potentialSites.push(combinedSite); }
            } catch (error) { console.error(`[Worker] Error combining pair (${i + 1}, ${j + 1}):`, error, r1, r2); }
        }
    }
    console.log(`[Worker] Generated ${potentialSites.length} potential combined sites.`);

    if (potentialSites.length === 0) {
        console.log("[Worker] --- Finished Site Combination (No pairs combined) ---");
        return [];
    }

    console.log("[Worker]  Clustering potential sites...");
    const finalSites = []; let unassignedSites = [...potentialSites];
    while (unassignedSites.length > 0) {
        let currentGroup = [unassignedSites.shift()]; let remainingSites = [];
        for (const site of unassignedSites) { if (currentGroup.some(member => areClose(site.x, member.x) && areClose(site.y, member.y) && areClose(site.z, member.z))) { currentGroup.push(site); } else { remainingSites.push(site); } }
        unassignedSites = remainingSites;
        let sumX = 0, sumY = 0, sumZ = 0; const refX = currentGroup[0].x, refY = currentGroup[0].y, refZ = currentGroup[0].z;
        for (const site of currentGroup) { sumX += adjustPeriodic(site.x, refX); sumY += adjustPeriodic(site.y, refY); sumZ += adjustPeriodic(site.z, refZ); }
        const avgSite = { x: sumX / currentGroup.length, y: sumY / currentGroup.length, z: sumZ / currentGroup.length };
        const norm = val => (((val % 1) + 1) % 1);
        finalSites.push({ x: norm(avgSite.x), y: norm(avgSite.y), z: norm(avgSite.z), count: currentGroup.length });
        console.log(`[Worker]   Cluster (Size ${currentGroup.length}): Avg=(${finalSites[finalSites.length - 1].x.toFixed(3)}, ${finalSites[finalSites.length - 1].y.toFixed(3)}, ${finalSites[finalSites.length - 1].z.toFixed(3)})`);
    }

    console.log(`[Worker] --- Finished Site Combination (${finalSites.length} sites) ---`);
    return finalSites;
}


// --- WORKER MESSAGE HANDLER ---
/**
 * Steps 2-4 of the pipeline: peak-finding, Harker analysis, site combination.
 *
 * mapResolution here is the grid the FFT actually produced, which is not
 * necessarily the one the user requested - every index computed below depends
 * on getting that right.
 */
function runAnalysisSteps(pattersonMap3D, crystalData, spaceGroups, mapResolution, harkerTolerance) {
    // Fallback for the case where the symmetry was resolved elsewhere and
    // calculatePattersonMap never ran in this worker.
    if (!lastExpansionSymmetry && crystalData?.spaceGroup) {
        try {
            const resolved = resolveSpaceGroupSetting(spaceGroups, crystalData.spaceGroup.number, crystalData.spaceGroup.name);
            const ops = resolved ? getExpansionOperators(resolved.setting) : null;
            if (resolved && ops) {
                lastExpansionSymmetry = {
                    ok: true, source: ops.source, opCount: ops.opCount,
                    settingSymbol: resolved.setting.symbol, matched: resolved.matched,
                    warnings: resolved.warnings.concat(ops.warnings),
                    description: ops.source === 'sym_ops'
                        ? `full operators (${ops.opCount}) of ${resolved.setting.symbol}`
                        : `rotations only (${ops.opCount}, no translations) of ${resolved.setting.symbol}`
                };
            }
        } catch (e) { /* reported by the map path already */ }
    }
    if (lastExpansionSymmetry?.warnings?.length) {
        postMessage({ type: 'symmetry_warning', payload: lastExpansionSymmetry });
    }

    postMessage({ type: 'status', payload: 'Finding peaks...' });
    const foundPeaks = findPeaks(pattersonMap3D, mapResolution, crystalData?.cell);

    postMessage({ type: 'status', payload: 'Analyzing Harker sections...' });
    const harkerAnalysisResults = analyzeHarkerPeaks(foundPeaks, crystalData, spaceGroups, mapResolution);

    postMessage({ type: 'status', payload: 'Consolidating sites...' });
    const consolidatedSites = combineSites(harkerAnalysisResults, harkerTolerance);

    let finalMessage = "Done.";
    const numCombined = consolidatedSites.length;
    const numPartial = harkerAnalysisResults.length;
    const numPeaks = foundPeaks.length;
    if (numCombined > 0) { finalMessage = `Done. Found ${numCombined} site(s).`; }
    else if (numPartial > 0) { finalMessage = `Done. Found ${numPartial} partial sites, but none combined.`; }
    else if (numPeaks > 0) { finalMessage = `Done. Found peaks, but no Harker matches.`; }
    else { finalMessage = `Done. No significant peaks found.`; }

    return { foundPeaks, harkerAnalysisResults, consolidatedSites, finalMessage };
}

self.onmessage = (e) => {
    const { type, payload } = e.data;
    lastExpansionSymmetry = null;
    
    if (type === 'CALCULATE') {
        try {
            const { crystalData, spaceGroups, mapResolution, harkerTolerance } = payload;

            // Step 1: Calculate the map by FFT. Fast enough on this single
            // worker thread that the old main-thread WebGPU path (which ran
            // the same O(res^3 * reflections) summation the CPU used to, just
            // spread over more cores) has been removed entirely.
            postMessage({ type: 'status', payload: `Transforming ${mapResolution}^3 map...` });
            const t0 = (typeof performance !== 'undefined' ? performance.now() : Date.now());
            const { map: pattersonMap3D, res: actualRes, dMin, sigma } =
                calculatePattersonMap(crystalData, spaceGroups, mapResolution);
            const t1 = (typeof performance !== 'undefined' ? performance.now() : Date.now());
            console.log(`[Worker] ${actualRes}^3 Patterson map by FFT in ${Math.round(t1 - t0)} ms ` +
                        `(dMin ${isFinite(dMin) ? dMin.toFixed(2) : '?'} A, peak sigma ${sigma.toFixed(2)} A).`);

            // Everything downstream must use the grid the FFT actually chose,
            // not the one that was requested.
            const { foundPeaks, harkerAnalysisResults, consolidatedSites, finalMessage } =
                runAnalysisSteps(pattersonMap3D, crystalData, spaceGroups, actualRes, harkerTolerance);

            postMessage({
                type: 'analysis_complete',
                payload: { pattersonMap3D, foundPeaks, harkerAnalysisResults, consolidatedSites, finalMessage,
                           mapResolution: actualRes, dMin, sigma,
                           symmetry: lastExpansionSymmetry }
            }, [pattersonMap3D.buffer]);   // transfer, don't clone: the map can be tens of MB

        } catch (error) {
            // Send errors back to the main thread
            console.error("[Worker] Pipeline Error:", error);
            postMessage({ type: 'error', payload: error.message || "An unknown worker error occurred." });
        }
    }

    else if (type === 'COMBINE_ONLY') {
        try {
            const { harkerAnalysisResults, harkerTolerance } = payload;

            postMessage({ type: 'status', payload: 'Re-consolidating sites...' });
            const consolidatedSites = combineSites(harkerAnalysisResults, harkerTolerance);

            let finalMessage = "Re-combine complete.";
             if (consolidatedSites.length > 0) { finalMessage = `Re-combine complete. Found ${consolidatedSites.length} site(s).`; }
             else { finalMessage = `Re-combine complete. No sites combined with this tolerance.`; }

            // --- Send back ONLY the updated consolidated sites ---
            // The main thread still *has* the map, peaks, and partial sites.
            postMessage({
                type: 'combine_complete', // Use a distinct type
                payload: {
                    consolidatedSites: consolidatedSites,
                    finalMessage: finalMessage
                }
            });

        } catch (error) {
             console.error("[Worker] Pipeline Error (COMBINE_ONLY):", error);
            postMessage({ type: 'error', payload: error.message || "An unknown worker error occurred during combine." });
        }
    }


};