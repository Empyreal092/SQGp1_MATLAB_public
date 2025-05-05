%standard 2D fourier transformation
function res=fft2_n(x)

NN = size(x);
if length(NN) == 2
    [Nx,Ny]=size(x);
    res=fft2(x)/(Nx*Ny); % correct normalization
else
    Nx = NN(end-1); Ny = NN(end);
    x = squeeze(x);
    res(1,1,:,:)=fft2(x)/(Nx*Ny); % correct normalization
end


end