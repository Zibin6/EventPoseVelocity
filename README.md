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

@ARTICLE{11559235,
  author={Liu, Zibin and Liang, Shunkun and Guan, Banglei and Shang, Yang and Yu, Qifeng and Zhao, Ji},
  journal={IEEE Transactions on Image Processing}, 
  title={A Geometric Framework for Absolute Pose and Velocity Estimation With Event Cameras}, 
  year={2026},
  volume={35},
  number={},
  pages={6714-6729},
  keywords={Cameras;Estimation;Vectors;Pose estimation;Optimization;Educational institutions;Event detection;Matrices;Timing;6-DOF;Event camera;pose estimation;velocity estimation;geometric solution},
  doi={10.1109/TIP.2026.3700937}}

