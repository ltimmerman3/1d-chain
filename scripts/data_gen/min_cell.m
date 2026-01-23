function cell_se = min_cell(Z_in, sigma_in, kappa, epsilon, cell_ref)

    boxes = linspace(1,20,39);
    energies = zeros(39,1);

    for i = 1:39
        lattice = boxes(i);
        dist = lattice / 2;
        S = msparc_1d_chain_paramOpt(dist, lattice, Z_in, sigma_in, kappa, epsilon, 2, 'GGA_PBE');
        energies(i) = S.Etotal;
    end

    xs = linspace(1,20,381);
    energy_curve = interp1(boxes, energies, xs, 'spline');

    % Find the minimum lattice constant on this continuous spline curve
    % Bounds are restricted to the data range to prevent wild extrapolation
    [~,opt_cell_idx] = min(energy_curve);
    opt_cell_spline = xs(opt_cell_idx);

    % Calculate squared error against the reference
    cell_se = (opt_cell_spline - cell_ref)^2;

    % plot(xs, energy_curve);

end