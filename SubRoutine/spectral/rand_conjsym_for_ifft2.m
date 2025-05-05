function [sym_M] = rand_conjsym_for_ifft2(n)
% make Gaussian random matrix that produce real ifft2 output.
% cf: https://www.mathworks.com/matlabcentral/answers/217750-make-ifft2-to-give-real-output

% Fill random examplary matrices
M = randn(n)+1i*randn(n);
% Apply symmetry to matrices
M(1,1)              = 0 ;
M(2:end,2:n/2)      = conj( rot90(M(2:end,n/2+2:end),2) );
M(1,2:n/2)          = conj( rot90(M(1,n/2+2:end),2) );
M(n/2+2:end,1)      = conj( rot90(M(2:n/2,1),2) );
M(n/2+2:end,n/2+1)  = conj( rot90(M(2:n/2,n/2+1),2) );
M(1,n/2+1)          = 2*randn;
M(n/2+1,n/2+1)      = 2*randn;
M(n/2+1,1)          = 2*randn;

sym_M = M;

end

