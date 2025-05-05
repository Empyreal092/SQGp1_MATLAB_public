function [rhs_k,u,v] = RHS_sqgp1(b_k,Ro,ikx_,iky_,K_,IK_,aa_filter)
%RHS_SQGP1 Summary of this function goes here
%   Detailed explanation goes here
phik = b_k.*IK_;
b_k_mean = b_k(1,1);
b_k(1,1)=0;

psix = ifft2_n(ikx_.*phik,'symmetric');
psiy = ifft2_n(iky_.*phik,'symmetric');
bx   = ifft2_n(ikx_.*b_k,'symmetric');
by   = ifft2_n(iky_.*b_k,'symmetric');

if Ro ~= 0
    b    = ifft2_n(b_k,'symmetric');
    bz   = ifft2_n(K_.*b_k,'symmetric');

    Rxk = K_.*aa_cut(fft2_n(b.*psix),aa_filter);
    Ryk = K_.*aa_cut(fft2_n(b.*psiy),aa_filter);
    bbz_k = aa_cut(fft2_n(b.*bz),aa_filter);
    Sk  = IK_.*bbz_k;
    
    u = -psiy.*(1+Ro*bz) + Ro*ifft2_n(iky_.*Sk + Ryk,'symmetric');
    v =  psix.*(1+Ro*bz) - Ro*ifft2_n(ikx_.*Sk + Rxk,'symmetric');
else
    u = -psiy;
    v =  psix;
end

rhs_k = -aa_cut(fft2_n( u.*bx + v.*by ),aa_filter);

end

