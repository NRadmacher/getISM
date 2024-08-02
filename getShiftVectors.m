function [shiftVector, save_name] = getShiftVectors(fname, dcname)

tmp_name        = strsplit(fname, '\');

if contains(tmp_name{end-1},".")
    date     = tmp_name{end-2};
else
    date     = tmp_name{end-1};
end
% date            = tmp_name{end-1};
img_name        = tmp_name{end};
img_name        = strsplit(img_name, '.');
img_name        = img_name{end-1};
img_name        = append(date,' ',img_name);

% [im_chan,~,im_posy,im_posx,~, head] = read_ISM(fname);
% 
% [head, ~, im_chan, im_posy, im_posx, ~] = PTU_ScanRead(fname, 0);
% 
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

%Dark count in count per second
[dc, ~] = get_DC(dcname,0);

%number of pixels of the detector
if(sum(head.HWInpChan_Enabled) > 23)
    n_pixl = 32;
    title_name = 'MPMT shift vectors';
    save_name = 'MPMT_shift_vectors';
    pix_time = head.ImgHdr_PixelTime;
    i = 9;
else
    n_pixl = 23;
    title_name = 'shift vectors';
    save_name = 'SPAD_shift_vectors';
%     im_chan = im_chan-9;
    pix_time = head.ImgHdr_PixelTime;
    i = 12;
end

if strcmp(head.CreatorSW_Name, 'SymPhoTime 64')
    n_pixl = 23;
    title_name = 'shift vectors';
    save_name = 'SPAD_shift_vectors_FlimBee';
    %timer per pixel in total
    pix_time = head.ImgHdr_MaxFrames * head.ImgHdr_TimePerPixel/1e3;
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
%        img = max(img - dc(pixl) * pix_time, 0);
       img = max(img -0, 0);
       
       imgs{pixl} = img;
end
disp("conversion to imges done")

shift_x_ic = zeros(n_pixl,1);
shift_y_ic = zeros(n_pixl,1);
% x is fist coordinat and y second 

if im_res < 0.03
    wd = 14;
else
    wd = 7;
end

center = imgs{i};
parfor j = 1:n_pixl
    
    %Phase correlarion to find shift between ism images(unshifted)
    %image i is the "center"
    
    [dx,dy] = image_corr(center,imgs{j}, wd);
%     [dx,dy] = phase_corr(center,imgs{j});
    shift_x_ic(j) =  dx;
    shift_y_ic(j) =  dy;
end

%% plots

xc  = -1.*(shift_x_ic);
yc  = -1.*(shift_y_ic);
%shift vector
sv_ic  = -[xc, yc];

h = figure;
ax = axes(h);
hold on
quiver(xc.*2, yc.*2 ,sv_ic(:,1),sv_ic(:,2),0, 'LineWidth', 2)  
numb = 0:n_pixl-1;
txt = string(numb);
plot(xc.*2, yc.*2,'xb','MarkerSize',10, 'LineWidth', 2, 'MarkerEdgeColor', 'green')
% plot(xc.*2, yc.*2 ,'o','MarkerSize',10, 'LineWidth', 2, 'MarkerEdgeColor', 'red')
set(ax,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1]);
ylabel('y shift [pixel]')
xlabel('x shift [pixel]')
XL = get(ax, 'XLim');
xl = XL(2) - XL(1);
textPosX = xc.*2+(0.025*xl);
textPosY =  yc.*2+(0.03*xl);
if n_pixl == 23
    textPosX(5) = -0.025*xl;
    textPosY(5) = -0.03*xl;
end
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
title(append(title_name,' image correlation'), 'FontSize', 15)
name = append(img_name, '_sv');
file_name = append(name,'.pdf');
exportgraphics(ax, file_name,'Resolution',600)

shiftVector = sv_ic.';
end

