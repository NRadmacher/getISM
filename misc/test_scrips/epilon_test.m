%% 
clear
close all
clc

%% Loadind data

fname       = 'W:\Niels\Messungen_Daten\240912_ism\lifetime.sptw\gattaBeadsRed_2.ptu';
caliname    = 'ismCallibration240912_ism.mat';

%nice colormap from TNT
cmap = cmap_greenFireBlue;

[im_time, ~, im_posx, im_posy, im_chan, head] = ScanRead(fname);

im_posx = double(im_posx);
im_posy = double(im_posy);

ind             = im_posy<=head.ImgHdr_PixY;
im_time         = im_time(ind);
im_chan         = im_chan(ind);
im_posx         = im_posx(ind);
im_posy         = im_posy(ind);
im_posx         = im_posx * head.ImgHdr_PixX;
im_time         = im_time./head.TTResult_SyncRate; % photon arrival in seconds
%% Parameters and magic numbers(please fix)

%number of pixels in recorded image
s_pixl_x    = head.ImgHdr_PixX;
s_pixl_y    = head.ImgHdr_PixY;

%scan resolution step per pixel in µm
IM_R = head.ImgHdr_PixResol;

%% seperate frames

% %select desierd frame binning
% binning = 20;
% binFrameStart = 73;
% binFrameEnd = 74;
% 
% %time frame  of imags
% sFrame  = binFrameStart*binning;
% eFrame  = binFrameEnd*binning;
% sTime   = sFrame * head.ImgHdr_FrameTime;
% etime   = eFrame * head.ImgHdr_FrameTime;
% ind     = im_time<etime & im_time>sTime;%9.5;
% 
% %only use photons in the selected time frame
% im_posx         = im_posx(ind);
% im_posy         = im_posy(ind);
% im_time         = im_time(ind);
% im_chan         = im_chan(ind);
%% ISM reasigment
fprintf('ISM reasigment ... ');

%shift vectors from file negativ sign is already included
calib       = open(caliname);
sv          = calib.shiftVector;
if calib.pixSize == IM_R
    psf = calib.psf;
else
    %calculate psf for current pixel size
    psf = calib.PSFfunc(calib.NA,calib.fd,calib.lamex,IM_R,calib.over);
end
shift_x     = sv(1, im_chan+1).';
shift_y     = sv(2, im_chan+1).';

%apply shift
ISM_posx    = im_posx + shift_x.*(calib.pixSize/IM_R);
ISM_posy    = im_posy + shift_y.*(calib.pixSize/IM_R);

clear shift_y shift_x;
fprintf('Done!\n');
%% plot

[optBinning] = getISMbinning(s_pixl_y,sv, 24);

ISM_binning = optBinning;
ISM_R = IM_R / ISM_binning;
%get image from photon stream
ismImg  = img_ps(ISM_posx,ISM_posy,s_pixl_x,s_pixl_y,ISM_binning);

f = figure;
ax = axes(f);
imagesc(ax, ismImg);
colormap(cmap)
set(ax,'DataAspectRatio', [1,1,1], ...
'PlotBoxAspectRatio',[1 1 1], ...
'XDir','normal', ...
'YDir','reverse');

%select Fov
rec = drawrectangle(ax);
cutPos = rec.Position;

%image with silder to adjust epsilon
handles.Image = imcrop(ismImg,cutPos);
handles.fig = f;
max_eps = 0.5;
min_eps = 0;
minor_step = (max_eps-min_eps)/100;
major_step = 10*minor_step;

%position and values of slider
handles.slider = uicontrol( 'Style','slider',...
                            'Position',[0 0 500 20],...
                            'SliderStep', [minor_step, major_step],...
                            'Min',min_eps,'Max',max_eps,'Value',0.07);

%action listener to update image after silder has changed
handles.Listener = addlistener(handles.slider,'Value','PostSet',...
    @(s,e) deconvolve(handles,ismImg,psf,cmap,cutPos));

%contrast button
handles.button = uicontrol('Style','pushbutton','String','Spot','Position',[505 0 40 20],'Callback',@(s,e) spot(handles));

%drwa initial image
colormap(cmap)
imagesc(ax,handles.Image);
set(gca,'DataAspectRatio', [1,1,1], ...
'PlotBoxAspectRatio',[1 1 1], ...
'XDir','normal', ...
'YDir','reverse');
colorbar(ax)


%% user functions
function deconvolve(handles, ISM_img,psf,cmap,cutPos)
    slider_value = get(handles.slider,'Value');
    image  = ISM_frw(ISM_img,psf,slider_value);
    image = imcrop(image,cutPos);
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

function spot(handles)
    sortedIntensities = sort(handles.Image(:));
    minIdx = round(0.75*numel(sortedIntensities));
    minval = sortedIntensities( minIdx );
    maxval = sortedIntensities( end );
    clim(handles.fig.CurrentAxes,[minval maxval])
end
