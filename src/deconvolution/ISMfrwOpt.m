function fISM = ISMfrwOpt(data, psf, eid, eps)
% ISMFRWOpt Performs optimal Fuoroer reweighting of an ISM image 
% date: pixel reassigned image
% psf experimental detection Point Spread Function of the microscope
% eid experimental Exitation Intensity Distribution of the microscope
% eps regulatization parameter
    if nargin < 3
        %regulation parameter to prevent divergence
        eps = 0.07;
    end
    %iamge size 
    [nx,ny] = size(data,[1 2]);
    if nx<1000 || ny<1000
        Nx = 1000;
        Ny = 1000;
    else
        Nx = nx;
        Ny = ny;
    end

    %offset to zero
    psf = psf(1:2:end,1:2:end);
    eid = eid(1:2:end,1:2:end);
    psf = psf - min(psf,[],'all');
    eid = eid - min(eid,[],'all');
    %remove possible NaNs
    psf(isnan(psf)) = 0.0;
    eid(isnan(eid)) = 0.0;
    %Fouier transform
    Fpsf = abs(fftshift(fft2(psf,Nx,Ny)));
    Feid = abs(fftshift(fft2(eid,Nx,Ny)));
    %normalise
    Fpsf = Fpsf./max(Fpsf,[],'all');
    Feid = Feid./max(Feid,[],'all');

    %Fourier transform ISM image
    Fimg = fftshift(fft2(data,Nx,Ny));
    fa = abs(Fimg); % modulus
    fp = angle(Fimg); % Phase angle
    %exoerimantal OTF
    OTFexp = Fpsf.*Feid;
    %calculate reweighting factor
    Ffull = Feid./(OTFexp + eps);
    Ffull(OTFexp < eps) = 0;
    %only reweight Fourier amplitude not phase
    fISM = abs((ifft2((Ffull.*fa.*exp(1i.*fp)))));
    fISM = fISM(1:nx,1:ny);
end