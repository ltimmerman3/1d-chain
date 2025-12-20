function Hx = evaluateExactExchangePotential(S,X,Hx)
spinor = 1;

if S.ACEFlag == 0
    ndrange = (1+(spinor-1)*S.N:spinor*S.N);

    Vexx = zeros(S.N,size(X,2));

    for i = 1:size(X,2)
        for j = 1:S.Nev
            psiq = S.psi_outer(ndrange,j);
            rhs = conj(psiq).*X(:,i);
            if S.exxmethod == 0 % Fourier space
                k_shift = [0, 0, 0];
                V_ji = poissonSolve_FFT(S,rhs,k_shift,S.const_by_alpha);
            else % Real space
                f = poisson_RHS(S,rhs,S.Screening_Op);
                V_ji = S.dA \ f;
            end
            Vexx(:,i) = Vexx(:,i) - S.occ_outer(j,1)*V_ji.*psiq;
        end
    end
    Hx = Hx + S.hyb_mixing*Vexx;
else
    ndrange = (1+(spinor-1)*S.N:spinor*S.N);
    Xi_times_psi = (transpose(S.Xi(ndrange,:))*X)*S.dx;
    Hx = Hx - S.hyb_mixing * (S.Xi(ndrange,:)*Xi_times_psi);
end

end