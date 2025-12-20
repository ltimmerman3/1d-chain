function [S] = scf(S)
fprintf('\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n');
fprintf(' Starting SCF iteration...\n');

% Electrostatic potential
S = poissonSolve(S);

% Exchange-correlation potential
S = exchangeCorrelationPotential(S);

% Effective potential
S = calculate_effective_potential(S);

% Initialize the mixing history vectors
S.X = zeros(S.N*S.nspden,S.MixingHistory);
S.F = zeros(S.N*S.nspden,S.MixingHistory);
S.mixing_hist_xkm1 = zeros(S.N*S.nspden,1);
S.mixing_hist_fkm1 = zeros(S.N*S.nspden,1);

% initialize history
S.mixing_hist_xkm1(1:S.N) = S.rho(:,1);
S.mixing_hist_xk = S.mixing_hist_xkm1;

if S.usefock == 1
    scf_tol_init = S.SCF_tol_init;
else
    scf_tol_init = S.SCF_tol;
end

S.lambda_f = 0.0;
S = scf_loop(S,scf_tol_init);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%          SCF LOOP Exact Exchange            %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exact exchange potential 
if S.usefock == 1
    S.usefock = S.usefock+1;
else
    return;
end

% Exchange-correlation potential
S = exchangeCorrelationPotential(S);

% Effective potential
S = calculate_effective_potential(S);

% Note: No need to update mixing_hist_xk. Done in mixing of first iteration.

% Exact exchange potential parameters
err_Exx = S.FOCK_TOL+1;
count_Exx = 1;

S.lambda_f = 0.0;
while(count_Exx <= S.MAXIT_FOCK)
    % Store orbitals and occupations for outer loop
    S.psi_outer = S.psi;                                                                                                                                                                                                                                                                                                 
    S.occ_outer = S.occ;
    
    if S.ACEFlag == 1
        S = ace_operator(S);
    end
    
    % Calculate estimation of Exact Exchange energy
    S = evaluateExactExchangeEnergy(S);
    
    Eexx_pre = S.Eex;

     % Start SCF with hybrid functional
    S = scf_loop(S,S.SCF_tol,count_Exx);
    
    % update exact exchange energy
    S.Etotal = S.Etotal + 2*S.Eex;
    S.Exc = S.Exc - S.Eex;
    % Calculate accurate Exact Exchange energy
    S = evaluateExactExchangeEnergy(S);
    % update exact exchange energy
    S.Etotal = S.Etotal - 2*S.Eex;
    S.Exc = S.Exc + S.Eex;
    
    err_Exx = abs(S.Eex - Eexx_pre)/S.n_atm;        
    fprintf(' Exx outer loop error: %.4e \n',err_Exx) ;
    if err_Exx < S.FOCK_TOL && count_Exx >= S.MINIT_FOCK
        break;
    end

    count_Exx = count_Exx + 1;
end
fprintf('\n Finished outer loop in %d steps!\n', (count_Exx - 1));
fprintf(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n');

% make sure next scf starts with normal scf
S.usefock = S.usefock+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = scf_loop(varargin)
S = varargin{1};
scf_tol = varargin{2};
count_Exx = -1;

if nargin == 3
	count_Exx = varargin{3};
elseif nargin > 3 || nargin < 2
	error('Too many input arguments.');
end

% SCF LOOP 
count_SCF = 1;

if count_Exx > 0
    S.MINIT_SCF = 1;
end

% time for one scf
t_SCF = 0; 

% start scf loop
while count_SCF <= S.MAXIT_SCF
    tic_eig = tic;

    rho_in = S.rho; % For mixing
    
    if S.usefock > 1
        fprintf(' Outer loop iteration number: %2d\n', count_Exx);
    end
    fprintf(' ========================= \n');
    fprintf(' Relaxation iteration: %2d \n SCF iteration number: %2d \n',S.Relax_iter,count_SCF);
    fprintf(' ========================= \n');

    [S.psi, S.EigVal] = eigSolver(S);

    % Solve for Fermi energy S.lambda_f and occupations
    S = occupations(S);

    % for density mixing, us rho_in to estimate total energy
    [S.Etotal,S.Eband,S.Exc,S.Exc_dc,S.Eelec_dc,S.Eent] = evaluateTotalEnergy(S);

    % Electron density
    S = electronDensity(S);

    % Error in SCF fixed-point iteration
    err = (norm(S.rho - rho_in))/(norm(S.rho));
    fprintf(' Error in SCF iteration: %.4e \n',err) ;

    % Mixing to accelerate SCF convergence
    % get rho_out
    rho_out = zeros(S.N*S.nspden,1);
    rho_out(1:S.N) = S.rho(:,1);
    % mixing
    [S, rho_out] = mixing(S,rho_out,S.mixing_hist_xk,count_SCF);
    % update
    S.rho(:,1) = rho_out(1:S.N);

    % update Veff
    S = poissonSolve(S);
    % Exchange-correlation potential
    S = exchangeCorrelationPotential(S);
    % Effective potential            
    S = calculate_effective_potential(S);

    scf_runtime = toc(tic_eig);
	t_SCF = t_SCF + scf_runtime; % add chebyshev filtering time to SCF time
    fprintf(' This SCF iteration took %.3f s.\n\n', scf_runtime);

    count_SCF = count_SCF + 1 ;

    if err < scf_tol && count_SCF > S.MINIT_SCF
        break;
    end
end

fprintf('\n Finished SCF iteration in %d steps!\n', (count_SCF - 1));
fprintf(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S] = calculate_effective_potential(S)
S.Veff = S.XCswitch*S.Vxc + S.phi;
S.Veff = real(S.Veff);
end