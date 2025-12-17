%% Generate data
clear;
clc

addpath('src')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameters
N = 41;
r1s = linspace(5, 13, N)';
r2s = linspace(5, 13, N)';

% Preallocate a 2D matrix for the energies
energies = zeros(N, N);

% Nested loops to calculate energies
for i = 1:N
    r1 = r1s(i);
    for j = 1:N
        r2 = r2s(j);
        
        % Run the simulation
        S = msparc_1d_chain_2dPES(r1, r2, 15, 0.01, 10, 'GGA_PBE');
        
        % Store result in the (i, j) position of the matrix
        energies(i, j) = S.Etotal; 
    end
end

% Prepare data for saving
% We use ndgrid to create coordinate matrices that match the loop order (i,j)
[R1_Grid, R2_Grid] = ndgrid(r1s, r2s);

% Create a 3-column array: [r1_value, r2_value, energy]
save_array = [R1_Grid(:), R2_Grid(:), energies(:)];

% Save to file
save('data/1d_32atoms_5:0.2:13r1r2_kappa0.01_epsilon10_GGA_PBE_energies.mat', 'save_array');

% Create Surface Plot
figure;
% Note: surf(X, Y, Z) assigns X to columns and Y to rows.
% Since energies(i,j) has r1 as rows and r2 as cols, we pass (r2s, r1s, energies).
surf(r2s, r1s, energies);
xlabel('r2 distance');
ylabel('r1 distance');
zlabel('Total Energy');
title('Potential Energy Surface');
colorbar;
view(3); % Set standard 3D view