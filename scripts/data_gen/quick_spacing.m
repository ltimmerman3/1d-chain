%% Generate data
clear;
clc
restoredefaultpath
addpath('../../src')

% S = msparc_1d_chain(10, 32, 0.01, 10, 'GGA_PBE');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ang_to_bohr = 1.88973;

boxes = linspace(1,30,59);
energies = zeros(59,1);
kappa = 0.01;
epsilon = 1.;
% Z_in = 3;
% sigma_in = 1.241;
% 
% for i = 1:39
%     lattice = boxes(i);
%     dist = lattice / 2;
%     S = msparc_1d_chain_paramOpt(dist, lattice, Z_in, sigma_in, kappa, epsilon, 2, 'GGA_PBE');
%     energies(i) = S.Etotal;
% end

for i = 1:59
    lattice = boxes(i);
    dist = lattice / 2;
    S = msparc_1d_chain(dist, -1, lattice, kappa, epsilon, 'GGA_PBE');
    energies(i) = S.Etotal;
end

xs = linspace(1,30,581);
energy_curve = interp1(boxes, energies, xs, 'spline');

%% Generate beautiful plots
% Define size in inches (e.g., single column width ~3.5 inches)
width = 7;
height = 4.5;
f = figure('Units', 'inches', 'Position', [1 1 width height]);
% f.OuterPosition  = [0,0,
% 18,14];
tiledlayout(1,1, 'Padding', 'compact', 'TileSpacing', 'compact');

hold on
% plot(dists, energies,'b-','LineWidth', 1.2);
plot(xs ./ ang_to_bohr, energy_curve, 'b-', 'LineWidth', 1.2);
% plot(dists2, energies2, 'r-', 'LineWidth', 1.2);

% xlabel('Interatomic distance (bohr)', 'Interpreter', 'latex', 'FontSize', 14);
% xlabel('R (bohr)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Cell (Ang)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$E$ (ha)', 'Interpreter', 'latex', 'FontSize', 14);

box on;               % Enclose the plot in a box
set(gca, ...
    'LineWidth', 1.5, ...       % Thicker axis lines
    'FontSize', 16, ...         % Tick label size
    'FontName', 'Helvetica', ...% Professional font (sans-serif is usually best for plots)
    'TickLabelInterpreter', 'latex'); % LaTeX for tick numbers
% set(gca, 'YTickLabel', num2str(energies(:,1), '%.4f'));

% xlim([min(boxes) max(boxes)]);
xlim([4 30 / ang_to_bohr]);
% xlim([4 30]);

% egend('Box 7', 'Box 23')

% exportgraphics(f, 'Sigma2.0.pdf', 'Resolution', 300);

%% Surface plot

% energies is (71*46) x 3: [Etotal, box, dist]
% boxes is 71x1 (or 1x71), dists is 46x1 (or 1x46)

% Etotal_vec = energies(:,1);
% 
% % Reshape back into a 71x46 grid matching (i,j)
% Egrid = reshape(Etotal_vec, [46, 71]).';   % transpose -> 71x46
% 
% % Create coordinate grids
% [DistGrid, BoxGrid] = meshgrid(dists, boxes);  % both 71x46
% 
% % Plot
% figure;
% surf(DistGrid, BoxGrid, Egrid, 'EdgeColor', 'none');
% xlabel('dists');
% ylabel('boxes');
% zlabel('E_{total}');
% title('E_{total}(box, dist)');
% colorbar;
% view(45, 30);

