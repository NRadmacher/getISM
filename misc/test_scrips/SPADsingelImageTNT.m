
fname ='W:\Niels\Messungen_Daten\250108_ism\aligment.sptw\GroupMeas_2\gattaBeadsRed_7_z100,60µm_1.ptu';

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

im_frame        = cumsum([1; diff(im_posy)<0]);

if isfield(head, 'ImgHdr_FrameNum')
    if max(im_frame)~= head.ImgHdr_FrameNum
        fprintf('Warning calculated numner of frames does not match header file!\n')
    end
end
nFrame = max(im_frame);
frameBinning = 1;
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
       im_f = im_frame(ind);
       
       img = img_ps(im_x, im_y, im_f, s_pixl_x, s_pixl_y, 1,nFrame,frameBinning);
       
       if ~isempty(img)
        imgs(:,:,pixl) = rescale(sum(img,3));
       end
end

TNTvisualizer(imgs, struct('title','conf','metadata',struct('pixelsize',0.05,'pixelsize_unit',[char(181) 'm'])));




