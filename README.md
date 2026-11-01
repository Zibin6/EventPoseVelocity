# A Geometric Framework for Absolute Pose and Velocity Estimation with Event Cameras

### Introduction

This repository provides code for the simultaneous estimation of absolute pose and velocity using event cameras. 
The package implements algorithms to recover the full 6-DoF motion state, specifically calculating:
Absolute Pose: Rotation and translation.
Velocity: Angular velocity and linear velocity.

Code contributions by: [Zibin Liu](https://github.com/Zibin6) and [Shunkun Liang](https://github.com/LiangSK98)

### func

- `AbsLin`			: Absolute Linear Solver.
- `AbsPol`			: Absolute Polynomial Solver.
- `VelLin`			: Linear Velocity Solver.
- `VelOpt`			: Optimized Velocity Solver.

### Quick Start

Run the script by typing `demo.m` in the MATLAB Command Window.

### Simulation Experiments

The test directory contains MATLAB scripts for the error analysis of absolute pose and velocity estimation under various conditions.

### Reference

```bibtex
@article{liu2026geometric,
  author  = {Liu, Zibin and Liang, Shunkun and Guan, Banglei and
             Shang, Yang and Yu, Qifeng and Zhao, Ji},
  title   = {A Geometric Framework for Absolute Pose and Velocity
             Estimation With Event Cameras},
  journal = {IEEE Transactions on Image Processing},
  volume  = {35},
  pages   = {6714--6729},
  year    = {2026},
  doi     = {10.1109/TIP.2026.3700937}
}

