function [V] = poissonSolve_FFT(S,rhs,k_shift,fft_const)
shift_ind = find(ismembertol(S.k_shift,k_shift,1e-8,'ByRows',true))+0;
if shift_ind < S.num_shift
    u = rhs .* S.neg_phase(:,shift_ind);
else
    u = rhs;
end
% u = reshape(u,S.Nx,S.Ny,S.Nz); % Not needed as dealing only with 1D
% const_by_alpha = zeros(S.Nx,S.Ny,S.Nz);
% const_by_alpha(:) = fft_const(shift_ind,:,:,:);

u_hat = fft(u);
const_by_alpha = zeros(S.N,1); % Just 1D
const_by_alpha(:) = fft_const(shift_ind,:); % 1D

V = ifftn(u_hat.*const_by_alpha);
if shift_ind < S.num_shift
    V = V(:) .* S.pos_phase(:,shift_ind);
else
    V = V(:);
end


V = real(V(:));

end