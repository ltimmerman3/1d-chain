energies = load('data/Z-2_Sigma-1/1d_32atoms_5:0.1:13dist_kappa0.01_epsilon10_GGA_PBE_energies.mat').save_array;

dists = energies(:,1);
energies = energies(:,2);

% Define size in inches (e.g., single column width ~3.5 inches)
width = 7;
height = 4.5;
f = figure('Units', 'inches', 'Position', [1 1 width height]);
% f.OuterPosition  = [0,0,18,14];
tiledlayout(1,1, 'Padding', 'compact', 'TileSpacing', 'compact');

plot(dists, energies,'b-','LineWidth', 1.2);

xlabel('Interatomic distance (bohr)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$E$ (ha)', 'Interpreter', 'latex', 'FontSize', 14);

box on;               % Enclose the plot in a box
set(gca, ...
    'LineWidth', 1.5, ...       % Thicker axis lines
    'FontSize', 16, ...         % Tick label size
    'FontName', 'Helvetica', ...% Professional font (sans-serif is usually best for plots)
    'TickLabelInterpreter', 'latex'); % LaTeX for tick numbers
set(gca, 'YTickLabel', num2str(energies, '%.4f'));

% xlim([min(dists) max(dists)]);
xlim([5 max(dists)]);

% exportgraphics(f, 'Sigma2.0.pdf', 'Resolution', 300);