function [FK] = isospectrum(Fkl)

% [FK] = ISOSPECTRUM(Fkl) 
%     Azimuthally integrate spectral field Fkl:  
%     FK = \int_0^{2\pi} Fkl K d\theta
%
%     Fkl fft2(F), with size(F) = [nx,nx], nx a power of 2, with no
%     dummy values removed and no shifting. 

[n,~] = size(Fkl);
kmax = (n)/2-1;

[kx_,ky_] = meshgrid([0:n/2 1-n/2:-1],[0:n/2 1-n/2:-1]);
K_ = sqrt(kx_.^2+ky_.^2);

FK = zeros([1 kmax+1]);

% Fkl = fftshift(Fkl);

for K = 0:kmax
    mask = (floor(K_+.5)-K)==0;
    FK(K+1) = sum(mask.*Fkl,'all');
end 

