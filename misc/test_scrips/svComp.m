%% Clear shit  up
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F)
clear
close all
clc

%%


fname ='W:\Niels\Messungen_Daten\240927_ism\lifetime.sptw\GroupMeas_2\gattaBeadsRed_18_z0,80µm_1.ptu';
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

%number of pixels of the detector
if isfield(head,'HWInpChan_Enabled')
    if(sum(head.HWInpChan_Enabled) > 23)
        n_pixl = 32;
        title_name = 'MPMT shift vectors';
        save_name = 'MPMT_shift_vectors';
        i = 9;
    else
        n_pixl = 23;
        title_name = 'shift vectors';
        save_name = 'SPAD_shift_vectors';
        i = 12;
    end
    pixShift = 0;
elseif isfield(head,'HW_InpChannels')
    n_pixl = head.HW_InpChannels;
    title_name = 'shift vectors';
    save_name = 's23_shift_vectors';
    pixShift = 1;
    i = 12;
    
end
if strcmp(head.CreatorSW_Name, 'SymPhoTime 64')
    n_pixl = 23;
    title_name = 'shift vectors';
    save_name = 'SPAD_shift_vectors_FlimBee';
    i = 12;
end
%% calculate shift vectors with image correlation and phase correlation
if min(im_chan) ==  0
    im_pix = im_chan+1;
else
    im_pix = im_chan;
end

im_posx = im_posx * s_pixl_x;
imgs = cell(n_pixl, 1);

for pixl=1:n_pixl
    %care THG
       ind  = (im_pix==pixl);
       
       im_x = im_posx(ind);
       im_y = im_posy(ind);
       
       img = img_ps(im_x, im_y, s_pixl_x, s_pixl_y, 1);
       
       if pixShift
           shift_img = circshift(img(1:2:end,:),-1,2);
           img(1:2:end,:) = shift_img;
       end
       
       imgs{pixl} = img;
end
disp("conversion to imges done")

%% SV
shiftXic = zeros(n_pixl,n_pixl);
shiftYic = zeros(n_pixl,n_pixl);

shiftXicS = zeros(n_pixl,1);
shiftYicS = zeros(n_pixl,1);
% x is fist coordinat and y second 

if im_res < 0.03
    wd = 40;
else
    wd = 20;
end

%just with center
center = imgs{i};
parfor j = 1:n_pixl
    %Phase correlarion to find shift between ism images(unshifted)
    %image i is the "center"
    [dx,dy] = image_corr(center,imgs{j}, wd);
    shiftXicS(j) =  dx;
    shiftYicS(j) =  dy;
end

%global shift
for i = 1:n_pixl
    for j = 1:i
        [dx,dy] = image_corr(imgs{i},imgs{j});
        shiftXic(i,j) =  dx;
        shiftXic(j,i) =  -dx;
        shiftYic(i,j) =  dy;
        shiftYic(j,i) =  -dy;
    end
end

%% plots center only

xcS  = shiftXicS;
ycS  = shiftYicS;
%shift vector
svIcS  = [xcS, ycS];

h = figure;
ax = axes(h);
hold on
quiver(xcS.*2, ycS.*2 ,-svIcS(:,1),-svIcS(:,2),0, 'LineWidth', 2)  
numb = 0:n_pixl-1;
txt = string(numb);
plot(xcS.*2, ycS.*2,'xb','MarkerSize',10, 'LineWidth', 2, 'MarkerEdgeColor', 'green')
% plot(xc.*2, yc.*2 ,'o','MarkerSize',10, 'LineWidth', 2, 'MarkerEdgeColor', 'red')
set(ax,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1]);
ylabel('y shift [pixel]')
xlabel('x shift [pixel]')
XL = get(ax, 'XLim');
xl = XL(2) - XL(1);
textPosX = xcS.*2+(0.025*xl);
textPosY =  ycS.*2+(0.03*xl);
text(textPosX, textPosY, txt)
grid on
box on
set(findall(h,'-property','FontSize'),'FontSize',12)
set(findall(h,'-property','LineWidth'),'LineWidth',1.5)
%make axis perfect square
limits = max([abs(ax.YLim), abs(ax.XLim)], [], 'all').*1.15;
ylim( [-limits, limits] );
xlim( [-limits, limits] );
set(ax,'YTick',get(ax,'XTick'));
title(append('center',' image correlation'), 'FontSize', 15)

%% gobal

xc  = mean(shiftXic);
yc  = mean(shiftYic);

%shift vector
svIc  = [xc; yc].';

h = figure;
ax = axes(h);
hold on
quiver(xc.*2, yc.*2 ,-svIc(:,1).',-svIc(:,2).',0, 'LineWidth', 2)  
numb = 0:n_pixl-1;
txt = string(numb);
plot(xc.*2, yc.*2,'xb','MarkerSize',10, 'LineWidth', 2, 'MarkerEdgeColor', 'green')
set(ax,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1]);
ylabel('y shift [pixel]')
xlabel('x shift [pixel]')
XL = get(ax, 'XLim');
xl = XL(2) - XL(1);
textPosX = xcS.*2+(0.025*xl);
textPosY =  ycS.*2+(0.03*xl);
text(textPosX, textPosY, txt)
grid on
box on
set(findall(h,'-property','FontSize'),'FontSize',12)
set(findall(h,'-property','LineWidth'),'LineWidth',1.5)
%make axis perfect square
limits = max([abs(ax.YLim), abs(ax.XLim)], [], 'all').*1.15;
ylim( [-limits, limits] );
xlim( [-limits, limits] );
set(ax,'YTick',get(ax,'XTick'));
title(append('Golbal',' image correlation'), 'FontSize', 15)
