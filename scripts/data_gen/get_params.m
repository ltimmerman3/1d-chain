%% Generate data
clear;
clc
restoredefaultpath
addpath('../../src')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ang_to_bohr = 1.88973;

cell_ref = 6.0 * ang_to_bohr;
Z_in = 2;

func = @(sigma) min_cell(Z_in, sigma, 0.01, 1, cell_ref);

lb = 0.1;
ub = 4.0;
[sol, fval] = fminbnd(func, lb, ub);

fprintf("Sigma: %f\nSquared Error: %f", sol, fval);