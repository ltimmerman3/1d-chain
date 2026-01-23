%% Generate data
clear;
clc

addpath('src')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start = 3;
finish = 13;
delta = finish - start;
spacing = 0.2;
N = floor(delta / spacing) + 1;
ang_to_bohr = 1.88973;

% Initialize arrays
energies = zeros(N,1);
boxes = linspace(start,finish,N) * ang_to_bohr;
dists = boxes ./ 2;
eig_vals = [];
% forces = zeros(N,2);

for i = 1:N
    cell = boxes(i);
    dist = dists(i);
    S = msparc_1d_chain(dist, -1, cell, 0.01, 10, 'GGA_PBE');
    energies(i) = S.Etotal;
    % forces(i,:) = S.force';
    eig_vals = [eig_vals; S.EigVal'];
end

%% Plotting

% Define size in inches (e.g., single column width ~3.5 inches)
width = 7;
height = 4.5;
f = figure('Units', 'inches', 'Position', [1 1 width height]);
% f.OuterPosition  = [0,0,18,14];
tiledlayout(1,1, 'Padding', 'compact', 'TileSpacing', 'compact');

plot(boxes, energies, 'b-', 'LineWidth', 1.2);

% xlabel('Interatomic distance (bohr)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Cell size (bohr)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$E$ (ha)', 'Interpreter', 'latex', 'FontSize', 14);

box on;               % Enclose the plot in a box
set(gca, ...
    'LineWidth', 1.5, ...       % Thicker axis lines
    'FontSize', 16, ...         % Tick label size
    'FontName', 'Helvetica', ...% Professional font (sans-serif is usually best for plots)
    'TickLabelInterpreter', 'latex'); % LaTeX for tick numbers
% set(gca, 'YTickLabel', num2str(energies(:,1), '%.4f'));

% xlim([min(boxes) max(boxes)]);
% xlim([min(boxes(I-5)) max(boxes)]);

% legend('Box 7', 'Box 23')

% exportgraphics(f, 'Sigma2.0.pdf', 'Resolution', 300);

% band structure at minimum energy
[M,I] = min(energies);
tmp_occ = [ones(1,2) zeros(1,10)];
occ_index = tmp_occ > 0.75;
occupied = eig_vals(I,occ_index);
occupied_is = linspace(1, length(occupied), length(occupied));
unoccupied = eig_vals(I,~occ_index);
unoccupied_is = linspace(occupied_is(end)+1, occupied_is(end)+1+length(unoccupied), length(unoccupied));

figure;
hold on;
plot(occupied_is, occupied, 'bo', 'MarkerFaceColor','none', 'LineStyle','none'); % blue circles
plot(unoccupied_is, unoccupied, 'r^', 'MarkerFaceColor','none', 'LineStyle','none'); % red triangles
hold off;

xlabel('Index (i)');
ylabel('$\epsilon_i$', 'Interpreter','latex');
axis square
grid on;
legend({'occupied','unoccupied'}, 'Location','northwest', 'FontSize', 14);