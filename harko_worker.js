// --- HEAVY CALCULATION FUNCTIONS (MOVED FROM MAIN SCRIPT) ---

/**
 * Solves a coordinate string like "u/2" based on a peak's (u,v,w).
 */
function solveCoordinate(solverString, peak) {
    if (solverString === '?') return '?';
    try {
        const sanitizedExpression = solverString.replace(/[^uvw\d\+\-\*\/\%\.\(\)\s]/g, '');
        const solverFunc = new Function('u', 'v', 'w', `return ${sanitizedExpression}`);
        const result = solverFunc(peak.u, peak.v, peak.w);
        if (typeof result !== 'number' || isNaN(result) || !isFinite(result)) { throw new Error(`Solver returned non-finite: ${result}`); }
        return (((result % 1) + 1) % 1).toFixed(3);
    } catch (e) { console.error(`Error solving: "${solverString}" for peak (${peak.u.toFixed(3)}, ${peak.v.toFixed(3)}, ${peak.w.toFixed(3)}):`, e); return 'err'; }
}

/**
 * Calculates the 3D Patterson map. This is the main blocking function.
 */
function calculatePattersonMap(crystalData, mapResolution) {
    try {
        const { cell, reflections } = crystalData;
        if (!reflections || reflections.length === 0) { throw new Error("No reflection data."); }
        if (!cell || !cell.a || !cell.b || !cell.c || isNaN(cell.a) || isNaN(cell.b) || isNaN(cell.c)) { throw new Error("Invalid cell data."); }

        const res = mapResolution;
        console.log(`[Worker] Calculating with resolution: ${res}`);

        const V = cell.a * cell.b * cell.c;
        if (!V || !isFinite(V) || V <= 0) { throw new Error(`Invalid volume: ${V}`); }

        const pattersonMap3D = new Float32Array(res * res * res);
        const PI2 = 2 * Math.PI;

        for (let iw = 0; iw < res; iw++) {
            for (let iv = 0; iv < res; iv++) {
                for (let iu = 0; iu < res; iu++) {
                    const u = iu / res, v = iv / res, w = iw / res;
                    let p = 0;
                    for (const r of reflections) {
                        if (isNaN(r.intensity) || isNaN(r.h) || isNaN(r.k) || isNaN(r.l)) continue;
                        p += r.intensity * Math.cos(PI2 * (r.h * u + r.k * v + r.l * w));
                    }
                    if (!isFinite(p)) { p = 0; }
                    pattersonMap3D[iw * res * res + iv * res + iu] = p / V;
                }
            }
        }
        
        console.log(`[Worker] Map calculated.`);
        return pattersonMap3D; // Return the result
    } catch (error) {
        console.error("[Worker] Map calc error:", error);
        throw error; // Re-throw to be caught by the main handler
    }
}

/**
 * Finds peaks in the calculated 3D map.
 */
function findPeaks(pattersonMap3D, mapResolution) {
    const res = mapResolution, map = pattersonMap3D;
    if (!map) { return []; }
    let peaks = [], maxVal = -Infinity, minVal = Infinity;
    for (let i = 0; i < map.length; i++) { const v = map[i]; if (isFinite(v)) { if (v > maxVal) maxVal = v; if (v < minVal) minVal = v; } }
    if (!isFinite(maxVal) || maxVal === minVal) { console.warn(`[Worker] Map flat/invalid. Skipping peaks.`); return []; }
    
    const threshold = minVal + (maxVal - minVal) * 0.15;
    
    for (let iw = 1; iw < res - 1; iw++) {
        for (let iv = 1; iv < res - 1; iv++) {
            for (let iu = 1; iu < res - 1; iu++) {
                const idx = iw * res * res + iv * res + iu, val = map[idx];
                if (val < threshold || !isFinite(val)) continue;
                let isMax = true;
                for (let dw = -1; dw <= 1 && isMax; dw++) {
                    for (let dv = -1; dv <= 1 && isMax; dv++) {
                        for (let du = -1; du <= 1 && isMax; du++) {
                            if (du === 0 && dv === 0 && dw === 0) continue;
                            const nv = map[(iw + dw) * res * res + (iv + dv) * res + (iu + du)];
                            if (isFinite(nv) && nv > val) { isMax = false; break; }
                        }
                    }
                }
                if (isMax) {
                    const normH = (maxVal > minVal) ? (val - minVal) / (maxVal - minVal) : 0;
                    peaks.push({ u: iu / res, v: iv / res, w: iw / res, height: normH });
                }
            }
        }
    }
    peaks.sort((a, b) => b.height - a.height);
    const foundPeaks = peaks.slice(0, 50);
    console.log(`[Worker] Found ${peaks.length} peaks. Kept ${foundPeaks.length}.`);
    return foundPeaks;
}

/**
 * Checks found peaks against Harker sections.
 */
function analyzeHarkerPeaks(foundPeaks, crystalData, spaceGroups, mapResolution) {
    let harkerAnalysisResults = [];
    if (!crystalData?.spaceGroup || foundPeaks.length === 0 || !spaceGroups) { console.log("[Worker] Skipping Harker."); return []; }
    const sgNumber = crystalData.spaceGroup.number; const sgData = spaceGroups[sgNumber];
    const gridSpacing = 1.0 / mapResolution; const tol = 1.5 * gridSpacing;
    console.log(`[Worker] Analyzing SG: ${sgNumber}. Tol: ${tol.toFixed(3)}`);
    if (!sgData?.harker_sections) { console.warn(`[Worker] No Harker data for SG ${sgNumber}.`); return []; }
    if (sgData.harker_sections.length === 0) { console.log(`[Worker] No Harker sections for SG ${sgNumber}.`); return []; }
    
    console.log(`[Worker] Found ${sgData.harker_sections.length} Harker sections.`);
    sgData.harker_sections.forEach((section, si) => {
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
self.onmessage = (e) => {
    const { type, payload } = e.data;
    
    // Only message type is 'CALCULATE'
    if (type === 'CALCULATE') {
        try {
            const { crystalData, spaceGroups, mapResolution, harkerTolerance } = payload;
            
            // --- Run the full pipeline, step by step ---
            
            // Step 1: Calculate Map
            postMessage({ type: 'status', payload: `Calculating ${mapResolution}^3 map...` });
            const pattersonMap3D = calculatePattersonMap(crystalData, mapResolution);
            
            // Step 2: Find Peaks
            postMessage({ type: 'status', payload: 'Finding peaks...' });
            const foundPeaks = findPeaks(pattersonMap3D, mapResolution);

            // Step 3: Analyze Harker
            postMessage({ type: 'status', payload: 'Analyzing Harker sections...' });
            const harkerAnalysisResults = analyzeHarkerPeaks(foundPeaks, crystalData, spaceGroups, mapResolution);

            // Step 4: Combine Sites
            postMessage({ type: 'status', payload: 'Consolidating sites...' });
            const consolidatedSites = combineSites(harkerAnalysisResults, harkerTolerance);

            // --- Determine final status message ---
            let finalMessage = "Done.";
            const numCombined = consolidatedSites.length; 
            const numPartial = harkerAnalysisResults.length; 
            const numPeaks = foundPeaks.length;
            if (numCombined > 0) { finalMessage = `Done. Found ${numCombined} site(s).`; }
            else if (numPartial > 0) { finalMessage = `Done. Found ${numPartial} partial sites, but none combined.`; }
            else if (numPeaks > 0) { finalMessage = `Done. Found peaks, but no Harker matches.`; }
            else { finalMessage = `Done. No significant peaks found.`; }

            // --- Send all results back at once ---
            postMessage({ 
                type: 'analysis_complete', 
                payload: {
                    pattersonMap3D: pattersonMap3D, // Send the map
                    foundPeaks: foundPeaks,
                    harkerAnalysisResults: harkerAnalysisResults,
                    consolidatedSites: consolidatedSites,
                    finalMessage: finalMessage
                }
            });

        } catch (error) {
            // Send errors back to the main thread
            console.error("[Worker] Pipeline Error:", error);
            postMessage({ type: 'error', payload: error.message || "An unknown worker error occurred." });
        }
    }
};