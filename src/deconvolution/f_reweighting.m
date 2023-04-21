function W_ISM_img = f_reweighting(data, IM_R)
    
    %double the sampling points to get reweighting at twice the k-vectors
    IM_R = IM_R/2;
    %iamge size 
    [nx,ny] = size(data);
    %maximal radial exten of image in µm
    r_max = max(nx,ny) / 2 * IM_R;
    %calculate excitation intensity distrubution and point spread function
    [eid, psf] = get_psf([0 ceil(sqrt(2)*r_max)], IM_R);
    [Nx,Ny] = size(psf);
    %cut eid and psf to image size/factor
    psf = psf(floor((Nx-nx))+(1:nx), floor((Ny-ny))+(1:ny));
    eid = eid(floor((Nx-nx))+(1:nx), floor((Ny-ny))+(1:ny));
    %Foureier transforamtion
    Fpsf = abs(fftshift(fft2(psf)));
    Feid = abs(fftshift(fft2(eid)));
    %normalize 
    Fpsf = Fpsf./max(Fpsf,[],'all');
    Feid = Feid./max(Feid,[],'all');
    
    Fprod = Fpsf.*Feid;
    %regulation parameter to prevent divergence 0.013
    eps = 0.05 * max(Fprod,[],'all');
    %Fourier transform ISM image
    Fimg = fftshift(fft2(data));
    fa = abs(Fimg); % modulus
    fp = angle(Fimg); % Phase angle
    %calculate reweighting factor
    Ffull = Fpsf./(eps + Fprod);
%     Ffull(Fprod < eps) = 0;
    %only reweight Fourier amplitude not phase
    W_ISM_img = abs(ifft2(ifftshift(Ffull.*fa.*exp(1i.*fp))));

end