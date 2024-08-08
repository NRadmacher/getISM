%% 
clear
% close all
clc

%% Loadind data

fname   = 'W:\Niels\Messungen_Daten\240807_ism\lifetime.sptw\gattaBeadsRed_6.ptu';
dcname  = 'W:\Niels\Messungen_Daten\240806_ism\lifetime.sptw\dc_WL_2.ptu';
caliname    = 'ismCallibration240807_ism.mat';

cmap =cmap_greenFireBlue;

[im_time, im_tcspc, im_posx, im_posy, im_chan, head] = ScanRead(fname);

im_posx = double(im_posx);
im_posy = double(im_posy);

ind             = im_posy<=head.ImgHdr_PixY;
im_tcspc        = im_tcspc(ind);
im_chan         = im_chan(ind);
im_posx         = im_posx(ind);
im_posy         = im_posy(ind);
im_time         = im_time(ind);
im_posx         = im_posx * head.ImgHdr_PixX;
im_time         = im_time./head.TTResult_SyncRate; % photon arrival in seconds

%Dark count in count per second
[dc, ~] = get_DC(dcname,0);

%% Parameters and magic numbers(please fix)

%number of pixels in recorded image
s_pixl_x    = head.ImgHdr_PixX;
s_pixl_y    = head.ImgHdr_PixY;
s_pixl      = s_pixl_x * s_pixl_y;
s_size      = [s_pixl_x s_pixl_y];

%scan resolution step per pixel in µm
IM_R = head.ImgHdr_PixResol;

%% ISM reasigment
disp("ISM reasigment")

%shift vectors from file negativ sign is already included
calib = open(caliname);
sv = calib.shiftVector;
shift_x     = sv(1, im_chan+1).';
shift_y     = sv(2, im_chan+1).';

ISM_posx    = im_posx + shift_x.*(0.05/IM_R);
ISM_posy    = im_posy + shift_y.*(0.05/IM_R);

%%

ISM_img  = img_ps(ISM_posx, ISM_posy, s_pixl_x, s_pixl_y,1);
W_ISM_img1 = ISM_frw(ISM_img,calib.psf);
W_ISM_img1 = W_ISM_img1(239:325,107:193);
handles.Image = ISM_img;
handles.fig = figure;
max_eps = 0.2;
min_eps = 0;
minor_step = (max_eps-min_eps)/100;
major_step = 10*minor_step;
handles.slider = uicontrol( 'Style','slider',...
                            'Position',[0 0 500 20],...
                            'SliderStep', [minor_step, major_step],...
                            'Min',min_eps,'Max',max_eps,'Value',0.07);
handles.Listener = addlistener(handles.slider,'Value','PostSet',@(s,e) deconvolve(handles, ISM_img, calib,cmap));
% axis off;
colormap(cmap)
imagesc(handles.Image);
% xlim([xsize *0.35 xsize*0.75])
% ylim([ysize *0.135 ysize*0.535])
set(gca,'DataAspectRatio', [1,1,1], ...
'PlotBoxAspectRatio',[1 1 1], ...
'XDir','normal', ...
'YDir','reverse');





function deconvolve(handles, ISM_img, calib,cmap)
    slider_value = get(handles.slider,'Value');
    image  = ISM_frw(ISM_img,calib.psf, slider_value);
    image = image(107:193,239:325);
    handles.Image = image;
    imagesc(handles.Image);
    set(gca,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1], ...
    'XDir','normal', ...
    'YDir','reverse');
    colormap(cmap)
    colorbar
    title(sprintf('slider value: %f',slider_value))
    drawnow
end

