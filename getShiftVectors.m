function [shiftVector, save_name] = getShiftVectors(fname, dcname)

%Image title and name from file name
tmp_name        = strsplit(fname, '\');
%find date in filepath
i = 1;
date = [];
datePat = digitsPattern(6);
while isempty(date)
    if i == size(tmp_name,2)
        ME = MException('Filemane:nodate','file %s \ncontains no date', fname);
        date = {'000000'};
        fprintf('Warning no date found in file name!\n')
        break;
    end
    date = extract(tmp_name{end-i},datePat);
    i = i+1;
end
date = date{1};

%check if file was recorded by Symphotime and add ws name to image name
if contains(tmp_name{end-i+2},'.')
    wsName  = strsplit(tmp_name{end-i+2},'.');
    date    = append(date,' ',wsName{1});
    %check if file is part of Group Measurment and add name of GM to image
    %name
    if (i-3) > 0
        gmName  = tmp_name{end-i+3};
        date    = append(date,' ',gmName);
    end
end

img_name        = tmp_name{end};
img_name        = strsplit(img_name, '.');
img_name        = img_name{end-1};
img_name        = append(date,' ',img_name);

[~, ~, im_posx, im_posy, im_chan, head] = ScanRead(fname);

im_posx     = double(im_posx);
im_posy     = double(im_posy);
im_frame    = cumsum([1; diff(im_posy)<0]);

if isfield(head, 'ImgHdr_FrameNum')
    if max(im_frame)~= head.ImgHdr_FrameNum
        fprintf('Warning calculated numner of frames does not match header file!\n')
    end
end
nFrame = max(im_frame);
frameBinning = 1;

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
if isfield(head,'HWInpChan_Enabled')
    if(sum(head.HWInpChan_Enabled) > 23)
        n_pixl = 32;
        title_name = 'MPMT shift vectors';
        save_name = 'MPMT_shift_vectors';
    else
        n_pixl = 23;
        title_name = 'shift vectors';
        save_name = 'SPAD_shift_vectors';
    end
    pixShift = 0;
elseif isfield(head,'HW_InpChannels')
    n_pixl = head.HW_InpChannels;
    title_name = 'shift vectors';
    save_name = 's23_shift_vectors';
    pixShift = 1;  
end
if strcmp(head.CreatorSW_Name, 'SymPhoTime 64')
    n_pixl = 23;
    title_name = 'shift vectors';
    save_name = 'SPAD_shift_vectors_FlimBee';
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
   im_f = im_frame(ind);
   
   img = img_ps(im_x, im_y, im_f, s_pixl_x, s_pixl_y, 1,nFrame,frameBinning);
   img = sum(img,3);
   if pixShift
       shift_img = circshift(img(1:2:end,:),-1,2);
       img(1:2:end,:) = shift_img;
   end
   
   imgs{pixl} = img;
end
disp("conversion to imges done")

% x is fist coordinat and y second 
shiftXic = zeros(n_pixl,n_pixl);
shiftYic = zeros(n_pixl,n_pixl);

if im_res < 0.03
    wd = 40;
else
    wd = 20;
end

%global shift
for i = 1:n_pixl
    for j = 1:i-1
        [dx,dy] = image_corr(imgs{i},imgs{j}, wd);
        shiftXic(i,j) =  dx;
        shiftXic(j,i) =  -dx;
        shiftYic(i,j) =  dy;
        shiftYic(j,i) =  -dy;
    end
end
xc  = mean(shiftXic);
yc  = mean(shiftYic);

%shift vector
svIc  = [xc; yc].';
%% plots

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
textPosX = xc.*2+(0.025*xl);
textPosY =  yc.*2+(0.03*xl);
% if n_pixl == 23
%     textPosX(5) = -0.025*xl;
%     textPosY(5) = -0.03*xl;
% end
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
file_name = append(name,'.png');
exportgraphics(ax, file_name,'Resolution',600)

shiftVector = svIc.';
end

