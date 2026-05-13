# Iterated Posterior Linearisation Filtering for Digital Carrier Synchronisation

This repository contains the MATLAB implementation of the algorithms proposed in the paper:

**"Iterated Posterior Linearisation Filtering for Digital Carrier Synchronisation"**  
by **Muyang Li** and **Ángel F. García-Fernández**  
Presented at the **2024 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI)**, held on **September 4, 2024**.

## Overview

This project focuses on filtering-based methods for digital carrier synchronisation, where carrier parameters such as phase and frequency are estimated from noisy nonlinear measurements.

The repository includes implementations of:

- Iterated Posterior Linearisation Filtering (IPLF)
- Probabilistic IPLF-related method
- Phase-Locked Loop (PLL)
- Sigma-Point Kalman Filtering (SKF)
- Unscented Kalman Filtering (UKF)

The code is organised by algorithm and is suitable for both reproducing the paper results and helping beginners learn Kalman filtering and nonlinear state estimation.
## Citation

If you find this repository helpful, please consider citing the following paper:

```bibtex
@inproceedings{li2024iterated,
  title     = {Iterated Posterior Linearisation Filtering for Digital Carrier Synchronisation},
  author    = {Li, Muyang and Garc{\'i}a-Fern{\'a}ndez, {\'A}ngel F.},
  booktitle = {Proceedings of the 2024 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI)},
  year      = {2024}
}
```
## Repository Structure

```text
.
├── IPLF/          # IPLF implementation and RMSE evaluation
├── IPLF_prob/     # Probabilistic IPLF-related implementation
├── PLL/           # PLL baseline method
├── SKF/           # Sigma-Point Kalman Filter implementation
├── UKF/           # Unscented Kalman Filter implementation
└── README.md
