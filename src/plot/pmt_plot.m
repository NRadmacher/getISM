function pmt_plot(im_chan, n_pixl, map)
%PMT_PLOT Plot count rate across MPMT pixels

chan_img   = histcounts(im_chan,n_pixl);
fiber_img   = zeros(64,1);
fiber_img(pix_code())   = chan_img;

xx = -9.05:2.3:7.05;
yy = -7.05:2.3:9.05;

[XX,YY] = meshgrid(xx, yy);

XX = repmat(XX(:), 1,4);
YY = repmat(YY(:), 1,4);

x_shift = 2.*[0, 0, 1, 1]; 
y_shift = 2.*[0, 1, 1, 0];

X = bsxfun(@plus, XX, x_shift);
Y = bsxfun(@plus, YY, y_shift);

numb = 1:64;
txt = string(numb);

figure
colormap(map)
patch(X.',Y.',fiber_img)
text(XX(:,1)+0.7, YY(:,1)+1, txt, 'Color','y')
set(gca,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1])
colorbar
title('MPMT')
grid("off")
axis("off")
end