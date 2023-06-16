function [W_ISM_img, psf] = f_reweighting(data, IM_R)
    %returns Fourier-reweithed ISM image and the PSF used for reweigthing

    %iamge size 
    [nx,ny] = size(data);
    %maximal radial exten of image in µm
    r_max = ceil(sqrt(2) * max(nx,ny) / 2 * IM_R);
    %calculate excitation intensity distrubution and point spread function
    [eid, psf] = get_psf([0 r_max], IM_R);
    [Nx,Ny] = size(psf);
    %cut eid and psf to image size/factor
    psf = psf(floor((Nx-nx)/2)+(1:nx), floor((Ny-ny)/2)+(1:ny));
    eid = eid(floor((Nx-nx)/2)+(1:nx), floor((Ny-ny)/2)+(1:ny));
    %Foureier transforamtion
    Fpsf = abs(fftshift(fft2(psf)));
    Feid = abs(fftshift(fft2(eid)));
    %normalize 
    Fpsf = Fpsf./max(Fpsf,[],'all');
    Feid = Feid./max(Feid,[],'all');

    %psf from file
%     psfData = matfile("PSF.mat");
%     psf = psfData.psf;
% %     [psf, ~,~] = FourierUpsampling(psf,2);
%     Fpsf = abs(fftshift(fft2(psf,nx,ny)));
%     Fpsf = Fpsf./max(Fpsf,[],'all');

    %cut off
%     [xx, yy] = meshgrid(1:nx, 1:ny);
%     u = zeros(size(xx));
%     km = 4*1/(0.64/(2*1.49));
%     u( (floor(xx-nx/2).^2 + floor(yy-ny/2).^2)<(km)^2 ) = 1;

%     Fpsf = Fpsf.*u;
%     Feid = Feid.*u;
%     
    Fprod = Fpsf.*Feid;
    %regulation parameter to prevent divergence 0.013. 0.03 for 1/psf
    eps = 0.03;
    %Fourier transform ISM image
    Fimg = fftshift(fft2(data));
    fa = abs(Fimg); % modulus
    fp = angle(Fimg); % Phase angle
    %calculate reweighting factor
    Ffull = Fpsf./(eps + Fprod);
%     Ffull = u./(Feid + eps* sqrt((floor(xx-nx/2).^2 + floor(yy-ny/2).^2))/km );
%     Ffull = 1./(Fpsf + eps);
    Ffull(Fpsf < eps) = 0;
    %only reweight Fourier amplitude not phase
    W_ISM_img = abs(ifft2(ifftshift(Ffull.*fa.*exp(1i.*fp))));

end