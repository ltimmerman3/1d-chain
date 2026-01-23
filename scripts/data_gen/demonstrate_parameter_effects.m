Z_in = 1;
n_atm = 2;

% Initialize arrays to loop over
sigmas = [0.5, 1.0, 2.0];
epsilons = [1 10 100];
kappas = [0.001 0.01 0.1 1];

% Initialize lattice and min_storage
lattices = linspace(1,20,39);
min_locs = zeros(10,2);

% Run for sigmas
for i = 1:numel(sigmas)
    energies = zeros(numel(lattices), 1);
    for j = 1:numel(lattices)
        epsilon = 10;
        kappa = 0.01;
        lattice = lattices(j);
        dist = lattice / n_atm;
        sigma_in = sigmas(i);
        S = msparc_1d_chain_paramOpt(dist, lattice, Z_in, sigma_in, kappa, epsilon, n_atm, 'GGA_PBE');
        energies(j) = S.Etotal;
    end
    xs = linspace(1,20,381);
    fine_energy_grid = interp1(lattices, energies, xs, 'spline');
    [low_energy, low_loc] = min(fine_energy_grid);
    min_locs(i,:) = [low_energy, xs(low_loc)];
end

% Run for epsilon
for i = 1:numel(epsilons)
    energies = zeros(numel(lattices), 1);
    for j = 1:numel(lattices)
        epsilon = epsilons(i);
        kappa = 0.01;
        lattice = lattices(j);
        dist = lattice / n_atm;
        sigma_in = 1.0;
        S = msparc_1d_chain_paramOpt(dist, lattice, Z_in, sigma_in, kappa, epsilon, n_atm, 'GGA_PBE');
        energies(j) = S.Etotal;
    end
    xs = linspace(1,20,381);
    fine_energy_grid = interp1(lattices, energies, xs, 'spline');
    [low_energy, low_loc] = min(fine_energy_grid);
    min_locs(3+i,:) = [low_energy, xs(low_loc)];
end

% Run for kappa
for i = 1:numel(kappas)
    energies = zeros(numel(lattices), 1);
    for j = 1:numel(lattices)
        epsilon = 10;
        kappa = kappas(i);
        lattice = lattices(j);
        dist = lattice / n_atm;
        sigma_in = 1.0;
        S = msparc_1d_chain_paramOpt(dist, lattice, Z_in, sigma_in, kappa, epsilon, n_atm, 'GGA_PBE');
        energies(j) = S.Etotal;
    end
    xs = linspace(1,20,381);
    fine_energy_grid = interp1(lattices, energies, xs, 'spline');
    [low_energy, low_loc] = min(fine_energy_grid);
    min_locs(6+i,:) = [low_energy, xs(low_loc)];
end

