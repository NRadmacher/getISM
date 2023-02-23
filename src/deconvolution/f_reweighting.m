function W_ISM_img = f_reweighting(data, IM_R)
    
    [nx,ny] = size(data);
    r_max = max(nx,ny) * IM_R;

    [eid, psf] = get_psf([0 ceil(sqrt(2)*r_max/4)], IM_R);
    [Nx,Ny] = size(psf);

    psf = psf(floor((Nx-nx/2)/2)+(1:nx/2), floor((Ny-ny/2)/2)+(1:ny/2));
    eid = eid(floor((Nx-nx/2)/2)+(1:nx/2), floor((Ny-ny/2)/2)+(1:ny/2));

    Fpsf = fftshift(fft2(psf));
    Feid = fftshift(fft2(eid));

    [Nx,Ny] = size(Fpsf);

    figure
    imagesc(abs(Fpsf));
    axis square
    
%     [X,Y]   = meshgrid( (1:Nx)-floor(Nx/2),(1:Ny)-floor(Ny/2));
%     [Xq,Yq] = meshgrid( (1:0.5:Nx)-floor(Nx/2),(1:0.5:Ny)-floor(Ny/2) );
%     Fpsf = interp2(X,Y,abs(Fpsf),Xq,Yq,'cubic');
%     Feid = interp2(X,Y,abs(Feid),Xq,Yq,'cubic');

    Fpsf = imresize(abs(Fpsf),2);
    Feid = imresize(abs(Feid),2);

    figure
    imagesc(abs(Fpsf));
    axis square

    Fprod = Fpsf.*Feid;

    eps = 0.001 * max(Fprod,[],'all');

    Fimg = fftshift(fft2(data));
    fa = abs(Fimg); % modulus
    fp = angle(Fimg); % Phase angle

    Fw = Fpsf./(eps+Fprod);
    Ffull = zeros(ISM_size);
    Fw_size = floor(size(Fw)./2);
    Fw_ind = [ISM_size(1)/2-Fw_size(1):ISM_size(1)/2+Fw_size(1);ISM_size(2)/2-Fw_size(2):ISM_size(2)/2+Fw_size(2)];
    Ffull(Fw_ind(1,:),Fw_ind(2,:)) = Fw;
    W_ISM_img = abs(ifft2(ifftshift(Ffull.*fa.*exp(1i.*fp))));

    N = 2^ceil(log2(max(nx,ny))+1); % size of image for Fourier transform
    % stich patches of the original image to form a periodic 'super-image'
    im = [data fliplr(data) data fliplr(data) data];
    im = [flipud(im); im; flipud(im); im; flipud(im); im; flipud(im)];
    % cut out the center part of the 'super-image' of size NxX
    [sx,sy] = size(im);
    im = im(floor((sx-N)/2)+(1:N),floor((sy-N)/2)+(1:N));
    % let the signal decay to zero towards the edges
    [X,Y] = meshgrid(1:N,1:N);
    X = (X-N/2-1);
    Y = (Y-N/2-1);
    R = sqrt(X.^2+Y.^2);
    W = 2-8*R/N; W = min(W,1); W = max(W,0);
    im = W.*im;
    % Do the Fourier transform
    fim = fftshift(fft2(im));
    fa = abs(fim); % modulus
    fp = angle(fim); % Phase angle
    % calculate k-vectors of the Fourier-transformed image
    fs = 1/0.05; % per µm
    k = fs/(N).*R;
    % kMax = 4 PI NA/lambda_em
    km = 32.3;
    ind = double(k<km); % Step function to clip everything outside km to zero
    % 'Ideal' MFT
    OTF = ind.*2/pi.*(acos(k./km)-k./km.*sqrt(1-(k./km).^2));
    % Re-weighting function
    ep = 0.015;
    w = ind./(OTF + ep.*k./km);
    % Back-transform the re-weighted FT
    tim = abs(ifft2(ifftshift(w.*fa.*exp(1i.*fp))));
    % cut out the original field of the data
    W_ISM_img = tim((N-nx)/2+(1:nx),(N-ny)/2+(1:ny));
end