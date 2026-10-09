# Flexible-Multibody-SSM

# Nonlinear Model Reduction of Flexible Multibody Systems via Spectral Submanifolds

[![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-blue.svg)](https://www.mathworks.com/)
[![SSMTool](https://img.shields.io/badge/SSMTool-v2.4-orange.svg)](https://github.com/haller-group/SSMTool-2.4)
[![COCO](https://img.shields.io/badge/COCO-Continuation-green.svg)](https://sourceforge.net/projects/cocotools/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the official MATLAB implementation for the paper: **"Model reduction of planar flexible multibody systems: spectral submanifold approach for floating frame of reference formulations"** submitted in the journal of Mechanism and Machine Theory.

## 📌 Overview

The dynamic analysis of highly flexible multibody systems is characterized by the complex nonlinear coupling between large rigid-body motions and continuous elastic deformations. While the **Floating Frame of Reference Formulation (FFRF)** accurately models these kinematics, it yields highly nonlinear equations of motion featuring deformation-dependent mass matrices and complex quadratic velocity forces.

Traditional linear projection methods (like CMS) fail to capture amplitude-dependent dynamics such as centrifugal stiffening and resonance shifts. This repository provides a **mathematically rigorous, nonlinear reduced-order modeling framework** based on the theory of **Spectral Submanifolds (SSMs)**.

### ✨ Core Scientific Contributions
1. **Implicit Tensor Extraction Algorithm:** We study an algorithm that completely bypasses the catastrophic symbolic inversion of the deformation-dependent mass matrix $\mathbf{M}(\mathbf{q})$. By applying implicit differentiation directly to the governing force balance at the equilibrium, the $\mathcal{O}(2)$ and $\mathcal{O}(3)$ Taylor tensors are extracted with absolute mathematical precision.
2. **Constrained DAE Reduction:** The framework successfully extends to closed-loop mechanisms (like the slider-crank) governed by index-3 Differential-Algebraic Equations (DAEs), utilizing a first-order descriptor state-space form.
3. **Transient Chaos via 4D SSMs:** We demonstrate the topological necessity of 4D manifolds for capturing highly complex, aperiodic energy exchanges (transient chaos) in lightly damped regimes.
4. **Nonlinear Visualization:** The code mathematically isolates the pure nonlinear manifold curvature, geometrically proving the exact capture of velocity-squared centrifugal stiffening forces.

---

## ⚙️ Installation & Dependencies

To run the codes in this repository, you need MATLAB (with the Symbolic Math Toolbox) and two external open-source packages:

### 1. SSMTool
This code relies on the **SSMTool** package developed by the Haller Group (ETH Zürich) for the computation of invariant manifolds.
*   **Download:** [SSMTool v2.4 (GitHub)](https://github.com/haller-group/SSMTool-2.4)
*   **Setup:** Clone or download the repository, extract it, and run `install.m` inside the SSMTool folder to add it to your MATLAB path.

### 2. COCO (Continuation Core)
For exact Full-Order Model (FOM) validation and numerical tracking of periodic orbits and saddle-node bifurcations, we use **COCO**.
*   **Download:** [COCO on SourceForge](https://sourceforge.net/projects/cocotools/)
*   **Setup:** Download and extract the package. In MATLAB, navigate to the COCO folder and run `startup.m`.

---

## 📂 Repository Structure

The repository is organized into three primary physical systems of increasing topological complexity. Each system includes a `build` script (for exact tensor extraction) and a `run` script (for SSM reduction, FRC plotting, and COCO validation).

### 1. Single Flexible Link (`/OneFlexibleLink_FFRF_SSM`)
*   `build_single_link_ffrf.m`: Extracts the exact tensors for the 4D state-space flexible pendulum.
*   `OneFlexibleLink_FFRF_SSM.m`: Computes the 2D SSM, extracts the softening backbone/FRC, and plots the pure nonlinear manifold.

### 2. Double Flexible Pendulum (`/Double_Pendulum`)
*   `build_double_link_ffrf.m`: Extracts exact tensors for the 8D state-space system.
*   `run_double_pendulum_2D.m`: Analyzes the moderate-energy regime (steady-state softening) using a 2D SSM.
*   `run_double_pendulum_4D_chaos.m`: Analyzes the lightly damped, high-energy regime. Constructs a 4D SSM to capture transient chaotic modal interactions and aperiodic beating.

### 3. Flexible Slider-Crank Mechanism (DAE System) (`/Slider_Crank`)
*   `build_slider_crank_dae.m`: Formulates the 7D descriptor state-space for the index-3 DAE system and extracts the constrained tensors.
*   `run_slider_crank_ssm.m`: Computes the constrained 2D SSM, perfectly capturing the severe saddle-node bifurcations and matching the index-1 Baumgarte-stabilized FOM continuations.

---

## 🚀 Quick Start / Usage

1. Open MATLAB and ensure both `SSMTool` and `COCO` are added to your path.
2. Navigate to the desired system folder (e.g., `Slider_Crank`).
3. Simply execute the run script:
   ```matlab
   run_slider_crank_ssm

    The script will automatically:

        Symbolically construct the FFRF model.

        Perform Exact Implicit Tensor Extraction (taking only a few seconds).

        Perform Linear Spectral Analysis.

        Extract the analytical Backbone and FRC via SSMTool.

        Run the exact numerical continuation via COCO.

        Generate highly academic, publication-ready plots (FRC comparisons, 3D Manifolds, and Pure Nonlinear Curvatures).

📖 References

If you use this code or the implicit tensor extraction algorithm in your research, please cite our paper:

    Draft Paper. [Journal Name], [Volume/Issue], [DOI/Link].

Foundational Literature:

    Haller, G., & Ponsioen, S. (2016). Nonlinear normal modes and spectral submanifolds: existence, uniqueness and use in model reduction. Nonlinear Dynamics, 86(3), 1493-1534.

    Jain, S., & Haller, G. (2022). How to compute invariant manifolds and their reduced dynamics in high-dimensional finite element models. Nonlinear Dynamics, 107(2), 1417-1450.

    Li, M., Jain, S., & Haller, G. (2023). Model reduction for constrained mechanical systems via spectral submanifolds. Nonlinear Dynamics, 111, 8881-8911.

✉️ Contact & Collaboration

This repository is maintained by Reza Nopour.
For questions regarding the mathematical formulation, the exact implicit tensor extraction algorithm, or potential collaborations, please feel free to reach out.

📧 Email: rnopour@aut.ac.ir
🏛 Institution: Department of Mechanical Engineering, Amirkabir University of Technology

Disclaimer: This code is provided as-is for academic and research purposes.
