function W_ISM_img = f_reweighting(data, IM_R)
    
    %iamge size 
    [nx,ny] = size(data);
    %maximal radial exten of image in µm
    r_max = max(nx,ny) / 2 * IM_R;
    factor = 1;
    %calculate excitation intensity distrubution and point spread function
    [eid, psf] = get_psf([0 ceil(sqrt(2)*r_max/factor)], IM_R);
    [Nx,Ny] = size(psf);
    %cut eid and psf to image size/factor
    psf = psf(floor((Nx-nx/factor)/factor)+(1:nx/factor), floor((Ny-ny/factor)/factor)+(1:ny/factor));
    eid = eid(floor((Nx-nx/factor)/factor)+(1:nx/factor), floor((Ny-ny/factor)/factor)+(1:ny/factor));
    %Foureier transforamtion
    Fpsf = fftshift(fft2(psf));
    Feid = fftshift(fft2(eid));
    %scale by factor
    Fpsf = imresize(abs(Fpsf),factor);
    Feid = imresize(abs(Feid),factor);
    %normalize 
    Fpsf = Fpsf./max(Fpsf,[],'all');
    Feid = Feid./max(Feid,[],'all');
    
    Fprod = Fpsf.*Feid;
    %regulation parameter to prevent divergence
    eps = 0.013 * max(Fprod,[],'all');
    %Fourier transform ISM image
    Fimg = fftshift(fft2(data));
    fa = abs(Fimg); % modulus
    fp = angle(Fimg); % Phase angle
    %calculate reweighting factor
    Ffull = Fpsf./(eps + Fprod);
%     Ffull(Fprod < eps) = 0;

    figure
    mesh(abs(Ffull));
    axis square
    %only reweight Fourier amplitude not phase
    W_ISM_img = abs(ifft2(ifftshift(Ffull.*fa.*exp(1i.*fp))));

end