% Clear workspace and command window
clear all; clc; close all;
addpath('func\')

% Set parameters for simulation
maxTime = 0.01;      % Maximum simulation time
eventNum = 2;        % Number of events per line
lineNum =6;          % Number of lines to simulate
samplNum=10;         % Number of sample
% Generate synthetic event data
[Line, evt, K, R_gt, t_gt, w_gt, v_gt] = generate_eventline(lineNum, samplNum, eventNum, maxTime);

% Perform absolute linear estimation
[R_lin, t_lin] = AbsLin(Line, evt(:, 1));

% Perform absolute polynomial estimation
[R_pol, t_pol] = AbsPol(Line, evt(:, 1), R_gt);

% Calculate linear velocity estimates
[w_lin, v_lin] = VelLin(Line, evt, R_gt, t_gt);

% Calculate optimized velocity estimates
w_init = [0;0;0];v_init = [0;0;0];
[w_opt, v_opt] = VelOpt(Line, evt, R_gt, t_gt, w_lin, v_lin);

% Rotation Matrices Table
fprintf('=== Rotation Matrix Estimation ===\n\n')
fprintf('%-15s | %-60s\n', 'Method', 'Rotation Matrix')
fprintf('---------------+------------------------------------------------------------\n')
fprintf('Ground Truth   | %s\n', formatMatrix(R_gt))
fprintf('AbsLin         | %s\n',  formatMatrix(R_lin))
fprintf('AbsPol         | %s\n\n', formatMatrix(R_pol))

% Translation Vectors Table
fprintf('=== Translation Vector Estimation ===\n\n')
fprintf('%-15s | %-25s\n', 'Method', 'Translation Vector')
fprintf('---------------+---------------------------\n')
fprintf('Ground Truth   | %s\n', formatVector(t_gt))
fprintf('AbsLin         | %s\n', formatVector(t_lin))
fprintf('AbsPol         | %s\n\n', formatVector(t_pol))

% Angular Velocities Table
fprintf('=== Angular Velocity Estimation ===\n\n')
fprintf('%-15s | %-25s\n', 'Method', 'Angular Velocity')
fprintf('---------------+---------------------------\n')
fprintf('Ground Truth   | %s\n', formatVector(w_gt))
fprintf('VelLin         | %s\n', formatVector(w_lin))
fprintf('VelOpt         | %s\n\n', formatVector(w_opt))

% Linear Velocities Table
fprintf('=== Linear Velocity Estimation ===\n\n')
fprintf('%-15s | %-25s\n', 'Method', 'Linear Velocity')
fprintf('---------------+---------------------------\n')
fprintf('Ground Truth   | %s\n', formatVector(v_gt))
fprintf('VelLin         | %s\n', formatVector(v_lin))
fprintf('VelOpt         | %s\n\n', formatVector(v_opt))

fprintf('================================================\n')
