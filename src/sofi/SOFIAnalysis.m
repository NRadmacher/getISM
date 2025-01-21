function [sof, im0, soffull] = SOFIAnalysis(im,ncum,ntime,nraster)

% soffull is a 5-dim array:
% 1st & 2nd dim = image cooridnates
% 3rd dim = correlation lag time 2^(n-1)
% 4th dim = cumulant order
% 5th dim = batch number

% im:       raw image
% ncum:     desierd cumulant order
% ntime:    width of window (in number for frames) which will be correlated in one batch
% nraster:  raster of redused pixel size

% max cumulant order:
if nargin<2 || isempty(ncum)
    ncum = 6;
end


% how often to coarsen time scale(number of lagtime doubelings):
if nargin<3 || isempty(ntime)
    ntime = ceil(log2(1e2/ncum));
else
    win = ntime;
    ntime = ceil(log2(ntime/ncum));
end
    
% how much to decrease pixe size
a = size(im,1);
b = size(im,2);
if nargin<4 || isempty(nraster)
    nraster = [0 0];
else
    nraster = [ceil(ceil((a-1)/2)*(nraster-1)) ceil(ceil((b-1)/2)*(nraster-1))];
    a = a + 2*nraster(1);
    b = b + 2*nraster(2);
end

if ncum>6
    error('SOFIAnalysis:argChk', 'Higher cumulants than 6th order are stupid')
end

if ischar(im)
    warning off
    im = double(fastTiff(im));
    warning on
end

%PLS FIX case 2^ntime*ncum > number of frames
if(size(im,3) ~= win )
    win = 2^ntime*ncum;
end
soffull = zeros(a,b,ntime,ncum-1,floor(size(im,3)/win));
im0 = zeros(a,b,floor(size(im,3)/win));
% tmp = zeros(a,b,win);
%going thought eatch batch of size win = ntime(old)
for k=1:floor(size(im,3)/win)
    if nargin>3 && ~(nraster(1)==0 &&nraster(2)==0)
        tmp = real(ifft2(ifftshift(padarray(fftshift(fft2(im(:,:,(k-1)*win+1:k*win))),nraster))));
    else
        tmp = im(:,:,(k-1)*win+1:k*win);
    end
    %for image mim
    im0(:,:,k) = sum(tmp,3);
    %compute deltaF = F - <F>_t
    tmp = tmp - repmat(mean(tmp,3),[1 1 size(tmp,3)]);
    %iterate over the number of time lage doubleings 
    for kk=1:ntime
        if kk>1
            tmp = (tmp(:,:,1:2:end-1)+tmp(:,:,2:2:end))/2;
        end
        %calc all cumulats up to ncum for given time lag kk
        for j=1:ncum-1
            soffull(:,:,kk,j,k) = cumulant0(tmp,j+1,1);
            disp([k kk j]);
        end
    end
    sof = squeeze(mean(sum(soffull,5),3));
%     tmp = sum(im0,3);
%     mim(cat(3,tmp,abs(sof)));
end

