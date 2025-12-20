function S = ace_operator(S)

S.Ns_occ = max(sum(S.occ_outer>1e-6));
S.Ns_occ = min(S.Ns_occ+S.EXXACEVal_state,S.Nev);
Ns = S.Ns_occ;

% Only gamma point code
S.Xi = zeros(S.nspinor*S.N,Ns);    % For storage of W and Xi

for spinor = 1:S.nspinor
    ndrange = (1+(spinor-1)*S.N:spinor*S.N);
    nsshift = 0; % No spin

    rhs = zeros(S.N,Ns);
    for i = 1:Ns
        rhs(:,i:Ns) = bsxfun(@times,S.psi_outer(ndrange,i:Ns),S.psi_outer(ndrange,i));
        V_i = zeros(S.N,Ns);
        for j = i:Ns
            if (S.occ_outer(i,1+nsshift) + S.occ_outer(j,1+nsshift) > 1e-4)
                if S.exxmethod == 0 % Fourier space
                    V_i(:,j) = poissonSolve_FFT(S,rhs(:,j),[0,0,0],S.const_by_alpha);
                else % Real space
                    f = poisson_RHS(S,rhs(:,j),S.Screening_Op);
                    V_i(:,j) = S.dA \ f;
                end
            end
        end
        S.Xi(ndrange,(i+1:Ns)) = S.Xi(ndrange,(i+1:Ns)) - S.occ_outer(i,1+nsshift)*bsxfun(@times,S.psi_outer(ndrange,i),V_i(:,(i+1:Ns)));
        S.Xi(ndrange,i) = S.Xi(ndrange,i) - bsxfun(@times,S.psi_outer(ndrange,(i:Ns)),V_i(:,(i:Ns))) * S.occ_outer((i:Ns),1+nsshift);
    end
    M = (transpose(S.psi_outer(ndrange,1:Ns))*S.Xi(ndrange,:))*S.dx;
    L = chol(-M);
    S.Xi(ndrange,:) = S.Xi(ndrange,:) / L; % Do it efficiently
end

end