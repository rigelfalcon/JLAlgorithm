

# Exponential Speedup of the Janashia-Lagvilava Matrix Spectral Factorization Algorithm

## 🚀 Project Overview

This project aims to achieve "exponential speedup" for high-dimensional matrix spectral factorization by non-commutative generalization of the core Janashia-Lagvilava algorithm. By extending the method to incorporate block matrix coefficients and a hierarchical processing strategy, this innovative approach significantly enhances computational efficiency for massive-dimensional matrices, commonly encountered in fields like neural data analysis, ultimately addressing a major bottleneck in traditional matrix spectral factorization.

**Paper Information:**
*   **Title:** Exponential Speedup of the Janashia-Lagvilava Matrix Spectral Factorization Algorithm
*   **ArXiv Preprint:** [arXiv:2503.02553v1](https://arxiv.org/abs/2503.02553)

## ✨ Core Innovations

*   **Non-Commutative Generalization:** Extends the central equation of the Janashia-Lagvilava method to the non-commutative case, allowing polynomial coefficients to be represented in block matrix form while preserving the equation's fundamental structure. Key functions involved include  `inverse_polynomial_block.m` (for block polynomial inversion).
*   **Hierarchical Block Processing:** Instead of incrementally making leading principal submatrices analytic step-by-step (e.g., 2×2, then 3×3), the new approach processes submatrices starting from the main diagonal and expanding in sizes that double sequentially (e.g., 2×2, 4×4, 8×8, 16×16, etc.).
*   **Real-time Processing Capability:** This method significantly boosts computational efficiency for large-scale matrices, extending the applicability to fields requiring real-time analysis of massive datasets.

## 🖥️ Current Code Status (In Development)

This repository contains the MATLAB implementation of the new Janashia-Lagvilava matrix spectral factorization algorithm.

*   **Functionality:** Successfully performs spectral factorization for matrices of dimension `d=1024`.
*   **Accuracy:** Achieves extremely high precision in the frequency domain, with an error of approximately `~3e-11`. The causality check in the time domain shows a residual error of approximately `~0.07`, which may be acceptable for many engineering applications.
*   **Performance:** In its purely serial execution on CPU, the `d=1024` case takes approximately `2` minutes. This is a substantial improvement over the `2.5` hours reported for the direct Janashia-Lagvilava algorithm in the paper.

## 🚧 To-Do & Future Optimizations

1.  **Optimize block hankel matrix inversion** 

2.  **GPU Acceleration (High Priority)**

3.  **Optimized CPU Parallelization Strategy** 

4.  **Code Modularization & Readability**

## ⚙️ How to Run

To set up the project and ensure all necessary paths are added to MATLAB, run the `setup.m` script:

```matlab
setup.m
```

Entrypoint for the main block-wise Janashia-Lagvilava algorithm:
`.\external\JLpackage\Janashia_Lagvilava_block.m`


## 🤝 Contributors & Motivation

This work is the result of the collaboration between **Professor Lasha Ephremidze** of Kutaisi International University, Georgia, and the team at the **Joint Cuba– China Laboratory for Frontier Research in Translational Neuroethology** (Ying Wang, Ronaldo García Reyes, and **Professor Pedro Valdes-Sosa**).

The primary motivation for this project is to address the significant computational challenges inherent in **large-scale brain causality inference**. By developing a highly efficient algorithm for accurate matrix spectral factorization, this research aims to provide a powerful tool for neuroscientists to analyze complex neural data and better understand the intricate dynamics of brain connectivity.

**Note:** This project is actively under development and optimization. Performance and results may improve with future code iterations.



