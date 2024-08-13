%% 
clear
% close all
clc

%% Loadind data
tic
% load('scimaps.mat');
%rgb values for color map black,blue,cyan,green,yellow,orange?,red,magenta
xl = 1:255/7:256;%like linespace but stepsize
yl = [[0 0 0]; [0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]; [1 0 1]];
lambda   = 1:256;
spectrum = interp1(xl, yl, lambda);
fname       = 'W:\Niels\Messungen_Daten\240801_ism\exAhi.sptw\otoferlin_Star635p_gold_ODp8_1.ptu';
dcname      = 'W:\Niels\Messungen_Daten\240806_ism\lifetime.sptw\dc_WL_2.ptu';
caliname    = 'ismCallibration240807_ism.mat';

tmp_name = strsplit(fname, '\');
date = tmp_name{end-1};
img_name = tmp_name{end};
img_name = strsplit(img_name, '.');
img_name = img_name{end-1};
img_name = append(date,' ',img_name);

%read ISM data from .ptu file
[im_chan,im_tcspc,im_posy,im_posx,im_time, head] = read_ISM(fname);
ind = im_posy<=head.ImgHdr_PixY;
im_chan = im_chan(ind);
im_posx = im_posx(ind);
im_posy = im_posy(ind);
im_posx = im_posx * head.ImgHdr_PixX;
im_pix = im_chan;
im_time = im_time./head.TTResult_SyncRate; % photon arrival in seconds

%remove time delay due to unsyncronised Multi Harps
% im_tcspc = remove_MHH_offset(im_tcspc,im_chan, head.max_bin, 90);

%Dark count in count per second
[dc, bin_dc] = get_DC(dcname,0);

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

% im_pix = im_pix - 9;

%shift vectors from file negativ sign is already included
%shift vectors from file negativ sign is already included
calib = open(caliname);
sv = calib.shiftVector;
n_pixl = size(sv,2);
shift_x     = sv(1, im_pix+1).';
shift_y     = sv(2, im_pix+1).';

ISM_posx    = im_posx + shift_x.*(0.05/IM_R);
ISM_posy    = im_posy + shift_y.*(0.05/IM_R);

%%

ISM_img  = img_ps(ISM_posx, ISM_posy, s_pixl_x, s_pixl_y,1);

handles.Image = ISM_img;
ysize = size(ISM_img,1);
xsize = size(ISM_img,2);
handles.fig = figure;
max_s = round(ysize*1.2);
min_s = round(ysize*0.8);
minor_step = 1/(max_s-min_s);
major_step = 10*minor_step;
handles.slider = uicontrol( 'Style','slider',...
                            'Position',[0 0 500 20],...
                            'SliderStep', [minor_step, major_step],...
                            'Min',min_s,'Max',max_s,'Value',ysize);
handles.Listener = addlistener(handles.slider,'Value','PostSet',@(s,e) deconvolve(handles, ISM_posx, ISM_posy, s_pixl_x, s_pixl_y));
% axis off;
colormap(hot)
imagesc(handles.Image);
% xlim([xsize *0.35 xsize*0.75])
% ylim([ysize *0.135 ysize*0.535])
set(gca,'DataAspectRatio', [1,1,1], ...
'PlotBoxAspectRatio',[1 1 1], ...
'XDir','normal', ...
'YDir','reverse');





function deconvolve(handles, ISM_posx, ISM_posy, s_pixl_x, s_pixl_y)
    slider_value = get(handles.slider,'Value');
    binning = slider_value/s_pixl_y;
    ISM_img  = img_ps(ISM_posx, ISM_posy, s_pixl_x, s_pixl_y,binning);
    ysize = size(ISM_img,1);
    xsize = size(ISM_img,2);
    handles.Image=ISM_img;
    imagesc(handles.Image);
%     xlim([xsize *0.35 xsize*0.75])
%     ylim([ysize *0.135 ysize*0.535])
    set(gca,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1], ...
    'XDir','normal', ...
    'YDir','reverse');
%     axis off;
    colormap(hot)
    title(sprintf('slider value: %f, binning: %f, y pixel: %f',slider_value, binning, ysize))
    drawnow
end

