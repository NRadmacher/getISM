function [W_ISM_img, psf] = f_reweighting(data, IM_R, over, Jpsf)
    %returns Fourier-reweithed ISM image and the PSF used for reweigthing
    
    %iamge size 
    [nx,ny] = size(data);

    if nx<1000 || ny<1000
        Nx = 700;
        Ny = 700;
    else
        Nx = nx;
        Ny = ny;
    end

    psf = Jpsf(1:2:end,1:2:end);
%     psf = Jpsf;
    psf = psf - min(psf,[],'all');
    psf(isnan(psf)) = 0.0;
    Fpsf = abs(fftshift(fft2(psf,Nx,Ny)));
    Fpsf = Fpsf./max(Fpsf,[],'all');
 
    %regulation parameter to prevent divergence 0.013. 0.03 for 1/psf
    eps = 0.1;
    %Fourier transform ISM image
    Fimg = fftshift(fft2(data,Nx,Ny));
    fa = abs(Fimg); % modulus
    fp = angle(Fimg); % Phase angle
    %calculate reweighting factor
    Ffull = 1./(Fpsf + eps);
    Ffull(Fpsf < eps) = 0;
    %only reweight Fourier amplitude not phase
    W_ISM_img = abs(ifftshift(ifft2(ifftshift(Ffull.*fa.*exp(1i.*fp)))));
    W_ISM_img = W_ISM_img(Nx/2:Nx/2+nx, Ny/2:Ny/2+ny);
end