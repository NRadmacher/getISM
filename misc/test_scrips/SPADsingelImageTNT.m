
fname ='W:\Niels\Messungen_Daten\240927_ism\lifetime.sptw\gattaBeadsRed_7.ptu';

[~, ~, im_posx, im_posy, im_chan, head] = ScanRead(fname);

im_posx = double(im_posx);
im_posy = double(im_posy);

%number of pixels in recorded image
s_pixl_x = head.ImgHdr_PixX;
s_pixl_y = head.ImgHdr_PixY;

ind = im_posy<=s_pixl_y;
im_chan = im_chan(ind);
im_posx = im_posx(ind);
im_posy = im_posy(ind);

%pixel size
im_res = head.ImgHdr_PixResol;

n_pixl = 23;

im_posx = im_posx * s_pixl_x;
imgs = zeros(s_pixl_x,s_pixl_y, n_pixl);
if min(im_chan) ==  0
    im_pix = im_chan+1;
else
    im_pix = im_chan;
end

for pixl=1:n_pixl
    %care THG
       ind  = (im_pix==pixl);
       
       im_x = im_posx(ind);
       im_y = im_posy(ind);
       
       img = img_ps(im_x, im_y, s_pixl_x, s_pixl_y, 1);
       
       if ~isempty(img)
        imgs(:,:,pixl) = img;
       end
end

TNTvisualizer(imgs, struct('title','conf','metadata',struct('pixelsize',0.05,'pixelsize_unit',[char(181) 'm'])));




