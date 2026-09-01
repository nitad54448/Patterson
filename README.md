# Laue Diagram Simulator

A web-based tool for simulating and visualizing Laue diffraction patterns for various crystal structures. This interactive simulator allows users to explore the relationship between crystal orientation, lattice parameters, and the resulting diffraction pattern in a Laue experiment.

The application is built entirely with HTML, CSS, and JavaScript, utilizing the `three.js` library for handling 3D calculations and rotations.

## Features

* **Crystal System Selection**: Supports all 7 crystal systems (Cubic, Tetragonal, Orthorhombic, Hexagonal, Trigonal, Monoclinic, Triclinic) and their corresponding Bravais lattices.
* **Adjustable Lattice Parameters**: Modify lattice constants (a, b, c) and angles (α, β, γ) to match specific materials. The inputs dynamically adapt to the selected crystal system.
* **Interactive Crystal Rotation**: Rotate the crystal in real-time along the X, Y, and Z axes using sliders or by clicking and dragging directly on the detector view.
* **Configurable X-Ray Source**: Set the minimum and maximum wavelength (λ) for the polychromatic X-ray beam used in the Laue method.
* **Flexible Detector Geometry**:
    * Switch between **back-scattering** and **transmission** setups.
    * Adjust the sample-to-detector distance.
    * Change the physical diameter of the detector.
* **Miller Indices Range**: Specify the range of (h, k, l) planes to be included in the simulation.
* **Spot Identification**: Click on any diffraction spot in the pattern to display its corresponding (h,k,l) Miller indices.
* **Responsive UI**: A resizable control panel allows for more or less space for the visualization.

## How to Use

1.  **Open the File**: Simply open the `Laue.html` file in any modern web browser. No server is needed.
2.  **Set Parameters**: Use the control panel on the left side to:
    * Choose a crystal system and Bravais lattice.
    * Enter the desired lattice parameters.
    * Set the crystal's initial orientation using the rotation sliders.
    * Define the wavelength range of the X-ray source.
    * Configure the detector's position, distance, and size.
3.  **Interact with the Pattern**:
    * The simulated Laue pattern will appear on the right.
    * **Click and drag** on the pattern to rotate the crystal interactively.
    * Use the **mouse wheel** while hovering over the pattern to zoom in or out (which adjusts the detector distance).
    * **Click** on a diffraction spot to see its (h,k,l) indices displayed in the info box.

## Technical Details

* **Frontend**: HTML, CSS, JavaScript
* **3D Graphics & Math**: `three.js` (v0.128.0) is used for vector and matrix operations to calculate the positions of reciprocal lattice points and handle crystal rotations.
* **Rendering**: The diffraction pattern is drawn on an HTML `<canvas>` element using the 2D rendering context.

## Attribution

* **Author**: NitaD, Université Paris-Saclay
* **Version**: 22.1 (Laue Sim)
* **Date**: August 2025

