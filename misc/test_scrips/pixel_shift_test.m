%% 
clear
close all
clc

%% Loadind data

fname       = 'W:\Florencia\Third Harmonic Project\SHG\brain_SG_004.ptu';

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

%% plot
%get image from photon stream
confImg  = img_ps(im_posx,im_posy,s_pixl_x,s_pixl_y,1);

f = figure;
ax = axes(f);
imagesc(ax, confImg);
colormap(cmap)
set(ax,'DataAspectRatio', [1,1,1], ...
'PlotBoxAspectRatio',[1 1 1], ...
'XDir','normal', ...
'YDir','reverse');

%select Fov
rec = drawrectangle(ax);
cutPos = rec.Position;

%image with silder to adjust epsilon
handles.Image = imcrop(confImg,cutPos);
handles.fig = f;
max_shift = 20;
min_shift = -20;
minor_step = 1/40;
major_step = 1/20;

%position and values of slider
handles.slider = uicontrol( 'Style','slider',...
                            'Position',[0 0 500 20],...
                            'SliderStep', [minor_step, major_step],...
                            'Min',min_shift,'Max',max_shift,'Value',0);

%action listener to update image after silder has changed
handles.Listener = addlistener(handles.slider,'Value','PostSet',...
    @(s,e) shift(handles,confImg,cmap,cutPos));

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
function shift(handles,img,cmap,cutPos)
    slider_value = get(handles.slider,'Value');
    fprintf('Shift set to %i\n',slider_value);
    shift_img = circshift(img(1:2:end,:),round(slider_value),2);
    img(1:2:end,:) = shift_img;
    image = imcrop(img,cutPos);
    imagesc(image);
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
