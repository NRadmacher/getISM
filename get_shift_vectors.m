function get_shift_vectors()

%terta 100 old
fname = 'D:\PHD\Data\2022\220309_terta_bead\tetra_beads_015.ptu';
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220302\tetra_beads_004.ptu';
dcname = 'D:\PHD\Data\2022\220106\cd_001.ptu';

%qdot 20nm
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220420\q_dot_003.ptu';
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220424\q_dot_008.ptu';

%good sv unit 08.07.2022
% fname = 'D:\PHD\Data\2022\220309\tetra_beads_015.ptu';

%good shift vektors with pinhole 220302 004
% fname   = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220314\homer_bassoon_006.ptu';

% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220726\qdots_em605nm_007.ptu';

%pc shift 
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\210310\Tubulin_Dylight488_002.ptu';

%sv for neurons sice 27.09.22 SPAD
fname = 'D:\PHD\Data\2022\220927\neurons_g1_cy2_syt1_one_005.ptu';
fname = 'D:\PHD\Data\2022\220513\tetra_beads_100nm_003.ptu';

%sc neurosn 01.01.23 PMT ?
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230112\neurons_g1_cy2_syt1_two_004.ptu';

%sv for neurons sice 28.09.22 200nm SPAD
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220928\neurons_g1_alexa488_psd95_one_002.ptu';

%for mpmt until 31.01.23
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\221117\neurons_g1_syt1_cy2_one_003.ptu';

%mpmt setup adjusted 01.02.23
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230202_tetra_beads\tetra_pmt_026.ptu';
%spad setup adjustmed 01.02.23
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230314_CMs_life\tetra_spad_017.ptu';

% dcname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230314_CMs_life\dc_spad_001.ptu'; %SPAD
% dcname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230203_g4_dc\dc_m15_pmt_001.ptu'; %PMT

%MITE Setup
fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230424_ISM_STORM\STORM_test_gatta_quant_beadsR_50nm_27p9_pix_3.ptu';
dcname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230418_ISM_STORM\beads_calibration_PSF_200mm_laser_off6.ptu';

tmp_name        = strsplit(fname, '\');
date            = tmp_name{end-1};
img_name        = tmp_name{end};
img_name        = strsplit(img_name, '.');
img_name        = img_name{end-1};
img_name        = append(date,' ',img_name);

[im_chan,~,im_posy,im_posx,~, head] = read_ISM(fname);

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
    save_name = 'MPMT_shift_vectors.m';
    pix_time = head.ImgHdr_PixelTime;
    i = 9;
else
    n_pixl = 23;
    title_name = 'shift vectors';
    save_name = 'SPAD_shift_vectors.m';
%     im_chan = im_chan-9;
    pix_time = head.ImgHdr_PixelTime;
    i = 12;
end

if strcmp(head.CreatorSW_Name, 'SymPhoTime 64')
    n_pixl = 23;
    title_name = 'shift vectors';
    save_name = 'SPAD_shift_vectors_mite.m';
    %timer per pixel in total
    pix_time = head.ImgHdr_MaxFrames * head.ImgHdr_TimePerPixel/1e3;
    i = 12;
end

%ISM scale factor 
alpha = 1;%(1 + 525/470);

%% calculate shift vectors with image correlation and phase correlation
im_pix = im_chan;

im_posx = im_posx * s_pixl_x;
imgs = cell(n_pixl, 1);

for pixl=1:n_pixl
       ind  = (im_pix==pixl-1);
       
       im_x = im_posx(ind);
       im_y = im_posy(ind);
       
       img = img_ps(im_x, im_y, s_pixl_x, s_pixl_y, 1);
       img = max(img - dc(pixl) * pix_time, 0);
       
       imgs{pixl} = img;
end
disp("conversion to imges done")

shift_x_pc = zeros(n_pixl,1);
shift_y_pc = zeros(n_pixl,1);

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
    shift_x_ic(j) =  dx;
    shift_y_ic(j) =  dy;
end

%% plots

xc  = -1.*(shift_x_ic)/alpha;
yc  = -1.*(shift_y_ic)/alpha;
% this is the shift vector for each channel no alpha needed here
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
text(xc.*2+(0.025*xl), yc.*2+(0.03*xl), txt)
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

inner = mean([norm(sv_ic(17,:)), norm(sv_ic(16,:)), norm(sv_ic(11,:)), norm(sv_ic(7,:)), norm(sv_ic(8,:)), norm(sv_ic(13,:))]);
disp(inner)
[pix_count,~] = histcounts(im_chan, n_pixl);
% hex_plot(pix_count, hot);
%% save 
shift_vector = sv_ic.';
save(save_name, 'shift_vector');
end

