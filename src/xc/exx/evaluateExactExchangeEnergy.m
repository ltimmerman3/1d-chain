function [S] = evaluateExactExchangeEnergy(S)
S.Eex = 0;
spinor = 1;

if S.ACEFlag == 0
    ndrange = (1+(spinor-1)*S.N:spinor*S.N);
    for i = 1:S.Nev
        for j = 1:S.Nev
            psiqi = S.psi_outer(ndrange,i);
            psikj = S.psi(ndrange,j);
            rhs = conj(psiqi) .* psikj;
            if S.exxmethod == 0 % Fourier space
                k_shift = [0, 0, 0];
                gij = poissonSolve_FFT(S,rhs,k_shift,S.const_by_alpha);
            else % Real space
                f = poisson_RHS(S,rhs,S.Screening_Op);
                gij = S.dA \ f;
            end

            S.Eex = S.Eex + S.occ_outer(i,1)*S.occ_outer(j,1)*real(sum(conj(rhs).*gij.*S.W));
        end
    end
else
    ndrange = (1+(spinor-1)*S.N:spinor*S.N);
    Psi_curr = S.psi(ndrange, 1:S.Ns_occ);
    Xi_curr  = S.Xi(ndrange, :);
    Overlaps = (Psi_curr' * Xi_curr) * S.dx;
    E_per_state = sum(abs(Overlaps).^2, 2);
    occ_vec = S.occ_outer(1:S.Ns_occ, 1);
    S.Eex = S.Eex + sum(occ_vec .* E_per_state);
end

S.Eex = -S.Eex*S.hyb_mixing/2;
fprintf(' Eex = %.8f\n', S.Eex);
end