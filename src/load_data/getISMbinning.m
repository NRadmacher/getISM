function [optBinning] = getISMbinning(n_pixl_y,shiftVectors,offPixel)
%GETISMBINNING Estimates the optimal pixel binning for ISM images to
%reduces scanning artifacts introduces by differnet detector geometries.
%Optimasiation is done by reducing the STD of an homgenious imput alonf the
%y direction (slow axis). 

%number do detector pixels
nPixl = size(shiftVectors,2);

%genrate homegenuois input
[xx, yy] = meshgrid(1:n_pixl_y,1:nPixl);
im_chan = reshape(yy,[nPixl*n_pixl_y, 1]);
im_posy = reshape(xx,[nPixl*n_pixl_y, 1]);
im_posx = ones(nPixl*n_pixl_y, 1);
im_frame = ones(nPixl*n_pixl_y, 1);

%some detectors have one pixel turend off
ind = im_chan ~= offPixel;
im_chan = im_chan(ind);
im_posx = im_posx(ind);
im_posy = im_posy(ind);
im_frame = im_frame(ind);

%apply shift vector for pixel reassigment
shift_y     = shiftVectors(2, im_chan);
im_posy = im_posy + shift_y.';

%
function ISMstd = ISMpixelShiftSTD(binning)
    [img, ~, ~] = img_ps(im_posx, im_posy, im_frame, 1, n_pixl_y,binning,1,1);
    % Pixels at the edge have less count due to ism pixel reassigment
    img = img(9:end-7);
    ISMstd = std(img);
end

%find minimum STD across all pixels
options = optimset('PlotFcns',@optimplotfval, 'TolX',1e-7);
optBinning = fminbnd(@ISMpixelShiftSTD,0.95,1.05, options);

% [img, ~, ~] = img_ps(im_posx, im_posy, 1, n_pixl_y,1);
% figure 
% plot(img(9:end-7))
% 
% [img, ~, ~] = img_ps(im_posx, im_posy, 1, n_pixl_y,optBinning);
% figure 
% plot(img(1:end))

end