%% Generate data
clear;
clc

addpath('src')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ang_to_bohr = 1.88973;
dists1 = linspace(0,16,81) *ang_to_bohr;
energies1 = zeros(81,1);
boxes = linspace(8)
cell = 8.0 * ang_to_bohr;

for i = 1:71
    cell = boxes(i);
    for j = 1:46
        dist = dists(j)*cell;
        S = msparc_1d_chain(dist, -1, cell, 0.01, 10, 'GGA_PBE');
        energies((i-1)*46+j,:) = [S.Etotal, boxes(i), dists(j)];
        forces((i-1)*46+j,:) = S.force';
    end
end


for i = 1:81
    dist = dists1(i);
    S = msparc_1d_chain(dist, -1, cell, 0.01, 10, 'GGA_PBE');
    energies1(i) = S.Etotal;
end