function get_shift_vectors()

% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220309\tetra_beads_015.ptu';
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220302\tetra_beads_004.ptu';
%qdot 20nm
fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220420\q_dot_003.ptu';
fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220424\q_dot_008.ptu';

%good sv unit 08.07.2022
fname = 'D:\PHD\Data\2022\220309\tetra_beads_015.ptu';

%good shift vektors with pinhole 220302 004
% fname   = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220314\homer_bassoon_006.ptu';

% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220726\qdots_em605nm_007.ptu';

%pc shift 
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\210310\Tubulin_Dylight488_002.ptu';

%sv for neurons sice 27.09.22 SPAD
fname = 'D:\PHD\Data\2022\220927\neurons_g1_cy2_syt1_one_005.ptu';

%sc neurosn 01.01.23 PMT ?
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230112\neurons_g1_cy2_syt1_two_004.ptu';

%sv for neurons sice 28.09.22 200nm SPAD
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220928\neurons_g1_alexa488_psd95_one_002.ptu';

%for mpmt until 31.01.23
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\221117\neurons_g1_syt1_cy2_one_003.ptu';

%mpmt setup adjusted 01.02.23
fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230202\tetra_pmt_026.ptu';

% dcname = 'D:\PHD\Data\2022\220106\cd_001.ptu'; %SPAD
dcname = 'D:\PHD\Data\2022\221117\dc_002.ptu'; %PMT

tmp_name        = strsplit(fname, '\');
date            = tmp_name{end-1};
img_name        = tmp_name{end};
img_name        = strsplit(img_name, '.');
img_name        = img_name{end-1};
img_name        = append(date,' ',img_name);

[im_chan,~,im_posy,im_posx,~, head] = read_ISM(fname);

%Dark count in count per second
[dc, ~] = get_DC(dcname,0);

%number of pixels of the detector
if(sum(head.HWInpChan_Enabled) > 23)
    n_pixl = 32;
    title_name = 'MPMT shift vectors';
    save_name = 'MPMT_shift_vectors.m';
else
    n_pixl = 23;
    title_name = 'SPAD shift vectors';
    save_name = 'SPAD_shift_vectors.m';
end

%number of pixels in recorded image
s_pixl_x = head.ImgHdr_PixX;
s_pixl_y = head.ImgHdr_PixY;

%ISM scale factor 
alpha = (1 + 525/470);

ind = im_posy<=s_pixl_y;
im_chan = im_chan(ind);
im_posx = im_posx(ind);
im_posy = im_posy(ind);


%% calculate shift vectors with image correlation and phase correlation
im_pix = im_chan;

im_posx = im_posx * s_pixl_x;
imgs = cell(n_pixl, 1);

for pixl=1:n_pixl
       ind  = (im_pix==pixl-1);
       
       im_x = im_posx(ind);
       im_y = im_posy(ind);
       
       img = img_ps(im_x, im_y, s_pixl_x, s_pixl_y, 1);
       img = max(img - dc(pixl) * head.ImgHdr_PixelTime, 0);
       
       J = wiener2(img,[20 20]);
       imgs{pixl} = J;
%        imgs{pixl} = img;
end
disp("conversion to imges done")

shift_x_pc = zeros(n_pixl,n_pixl);
shift_y_pc = zeros(n_pixl,n_pixl);

shift_x_ic = zeros(n_pixl,n_pixl);
shift_y_ic = zeros(n_pixl,n_pixl);

% x is fist coordinat and y second 
for i = 1:n_pixl
    for j = 1:n_pixl
        
        %Phase correlarion to find shift between ism images(unshifted)
        %image i is the "center"
        [dx,dy] = phase_corr(imgs{i},imgs{j});
        shift_x_pc(i,j) =  dx;
        shift_x_pc(j,i) =  -dx;
        shift_y_pc(i,j) =  dy;
        shift_y_pc(j,i) =  -dy;
        
        [dx,dy] = image_corr(imgs{i},imgs{j});
        shift_x_ic(i,j) =  dx;
        shift_x_ic(j,i) =  -dx;
        shift_y_ic(i,j) =  dy;
        shift_y_ic(j,i) =  -dy;
    end
end
xc  = -1.*mean(shift_x_pc);
yc  = -1.*mean(shift_y_pc);
% this is the shift vector for each channel
sv_pc  = -[xc; yc]./alpha;

%% plots
figure
axis equal
hold on
quiver(xc, yc ,sv_pc(1,:),sv_pc(2,:),0, 'LineWidth', 2)  
numb = 0:n_pixl-1;
txt = string(numb);
plot(xc, yc,'xb','MarkerSize',10, 'LineWidth', 2 ,'MarkerEdgeColor', 'green')
text(xc+0.2, yc, txt)
set(gca,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1]);
ylabel('y shift [pixel]')
xlabel('x shift [pixel]')
grid on
box on
title(append(title_name,' pc'))

xc  = -1.*mean(shift_x_ic);
yc  = -1.*mean(shift_y_ic);
% this is the shift vector for each channel
sv_ic  = -[xc; yc]./alpha;

h = figure;
ax = axes(h);
axis equal
hold on
quiver(xc, yc ,sv_ic(1,:),sv_ic(2,:),0, 'LineWidth', 2)  
numb = 0:n_pixl-1;
txt = string(numb);
plot(xc, yc,'xb','MarkerSize',10, 'LineWidth', 2, 'MarkerEdgeColor', 'green')
text(xc+0.2, yc, txt)
set(ax,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1]);
ylabel('y shift [pixel]')
xlabel('x shift [pixel]')
grid on
box on
title(append(title_name,' pc'))
name = append(img_name, '_sv');
file_name = append(name,'.png');
exportgraphics(ax, file_name,'Resolution',600)

% pixel_img   = histcounts(im_chan,0:n_pixl);
% hex_plot(pixel_img,hot)

%% save 
shift_vector = sv_ic;
save(save_name, 'shift_vector');
end

