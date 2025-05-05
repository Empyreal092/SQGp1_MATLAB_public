%standard 2D fourier inverse transformation
function res=ifft2_n(x,sym_flag)

arguments
    x
    sym_flag='notsym';
end

NN = size(x);
if length(NN) == 2
    [Nx,Ny]=size(x);
    if strcmp(sym_flag,'symmetric')
        res=ifft2(x,'symmetric')*(Nx*Ny); % correct normalization
    else
        res=ifft2(x)*(Nx*Ny); % correct normalization
    end
else
    Nx = NN(end-1); Ny = NN(end);
    x = squeeze(x);
    if strcmp(sym_flag,'symmetric')
        res=ifft2(x,'symmetric')*(Nx*Ny); % correct normalization
    else
        res=ifft2(x)*(Nx*Ny); % correct normalization
    end
end



end