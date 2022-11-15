function [img, lin_ind, im_size] = img_ps(x, y, x_pixel, y_pixel, scale_factor)
%IMG_PS Calculate 2D image form x,y photon stream data with subpixel scale
%factor by accumulating pixels
%The amount of pixels per dimension is i_pixel * scale_factor
% img       Intensity image
% lin_ind   linear index of evey recordet photon
% im_size   Size of Image after scaling, x,y

x_size      = round(x_pixel * scale_factor);
y_size      = round(y_pixel * scale_factor);
im_size     = [y_size x_size];
img_ind_x   = discretize(x,x_size);
img_ind_y   = discretize(y,y_size);

lin_ind     = sub2ind(im_size, img_ind_y, img_ind_x);

%order of y,x to fit the convention of imagesc()
img_ind     = [img_ind_y img_ind_x];
ind         = all(~isnan(img_ind),2);
img_ind     = img_ind(ind,:);

img         = accumarray(img_ind, 1);
end

