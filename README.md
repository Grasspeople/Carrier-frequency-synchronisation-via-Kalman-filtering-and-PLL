# Iterated Posterior Linearisation Filtering for Digital Carrier Synchronisation

This repository contains the official implementation of the algorithms proposed in the paper:

**"Iterated Posterior Linearisation Filtering for Digital Carrier Synchronisation"**  
by **Muyang Li** and **Ángel F. García-Fernández**  
Presented at the **2024 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI)**, held on **September 4, 2024**.

## Overview

Digital carrier synchronisation is a fundamental problem in communication and signal processing systems. It aims to estimate carrier-related parameters, such as phase and frequency, from noisy nonlinear measurements.

This repository provides MATLAB implementations of several filtering-based algorithms for digital carrier synchronisation, including:

- Iterated Posterior Linearisation Filtering (IPLF)
- Probabilistic IPLF-related implementation
- Phase-Locked Loop (PLL)
- Sigma-Point Kalman Filtering (SKF)
- Unscented Kalman Filtering (UKF)

The code is organised in a simple and modular way, making it suitable not only for reproducing the results of the paper, but also for beginners who want to learn and compare different Kalman filtering algorithms.

## Suitable for Beginners

This repository is also suitable for students and researchers who are new to Kalman filtering and nonlinear state estimation.

The implementation includes several commonly used filtering approaches, allowing beginners to:

- Understand the basic workflow of nonlinear filtering algorithms
- Compare IPLF, SKF, UKF, and PLL-based methods
- Learn how sigma points are generated and used
- Study how noisy measurements and ground-truth signals are generated
- Observe how filtering performance is evaluated using RMSE
- Use MATLAB scripts as simple examples for learning and modification

The code structure is relatively clear and each algorithm is placed in a separate folder, which makes it easier to follow and experiment with.

## Repository Structure

```text
.
├── IPLF/
│   ├── IPLF.m
│   ├── IPLF_RMSE.m
│   ├── draw_IPLF.m
│   ├── generate_sigma_point_IPLF.m
│   ├── generate_truth_IPLF.m
│   └── prob_drawing.m
│
├── IPLF_prob/
│   ├── IPLF_pro.m
│   ├── generate_sigma_point_IPLF_prob.m
│   ├── generate_truth_IPLF_prob.m
│   └── prob_drawing.m
│
├── PLL/
│   ├── PLL_RMSE.m
│   └── generate_truth_PLL.m
│
├── SKF/
│   ├── SKF_demo/
│   ├── PLL_RMSE.m
│   ├── SKF.m
│   ├── SKFfiltering_RMSE.m
│   ├── draw_filtered.m
│   └── generate_truth_measurement.m
│
├── UKF/
│   ├── UKF.m
│   ├── UKF_1.m
│   ├── UKF_RMSE.m
│   ├── calculate.m
│   ├── calculate2.m
│   ├── draw_UKF.m
│   ├── generate_sigma_point.m
│   └── generate_truth_UKF.m
│
└── README.md

}
## Citation

If you find this repository helpful, please consider citing the following paper:

```bibtex
@inproceedings{li2024iterated,
  title     = {Iterated Posterior Linearisation Filtering for Digital Carrier Synchronisation},
  author    = {Li, Muyang and Garc{\'i}a-Fern{\'a}ndez, {\'A}ngel F.},
  booktitle = {Proceedings of the 2024 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI)},
  year      = {2024}
