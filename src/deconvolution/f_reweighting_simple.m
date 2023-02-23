function W_ISM_img = f_reweighting_simple(data)
%F_REWEIGHTING_SIMPLE calclulate fourier reweighted ISM image modifiey
%after Gregor, I., Spiecker, M., Petrovsky, R. et al. Rapid nonlinear 
%image scanning microscopy. Nat Methods 14, 1087–1089 (2017).
    
    %iamge size 
    [nx,ny] = size(data);
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
    figure
    mesh(w);
    axis square
    % Back-transform the re-weighted FT
    tim = abs(ifft2(ifftshift(w.*fa.*exp(1i.*fp))));
    % cut out the original field of the data
    W_ISM_img = tim((N-nx)/2+(1:nx),(N-ny)/2+(1:ny));
end