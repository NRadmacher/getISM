function get_shift_vectors()

% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220309\tetra_beads_015.ptu';
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220302\tetra_beads_004.ptu';
fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220420\q_dot_003.ptu'; %qdot 20nm
fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220424\q_dot_008.ptu';

%good sv unit 08.07.2022
fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220309\tetra_beads_015.ptu';

% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220726\qdots_em605nm_007.ptu';

%pc shift 
fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\210310\Tubulin_Dylight488_002.ptu';
%sv for neurons sice 27.09.22
fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220927\neurons_g1_cy2_syt1_one_005.ptu';
dcname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220106\cd_001.ptu';
%good shift vektors with pinhole 220302 004
% fname   = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220314\homer_bassoon_006.ptu';

tmp_name        = strsplit(fname, '\');
date            = tmp_name{end-1};
img_name        = tmp_name{end};
img_name        = strsplit(img_name, '.');
img_name        = img_name{end-1};
img_name        = append(date,' ',img_name);

[im_chan,~,im_posy,im_posx,~, head] = read_ISM(fname);

det_y =-1.*[-2,-1,0,1,2,-3/2,-1/2,1/2,3/2,-2,-1,0,1,2,-3/2,-1/2,1/2,3/2,-2,-1,0,1,2];
A = ones(1,5).*sqrt(3);
B = zeros(1,5);
det_x = 1.*[A,A(1:end-1)./2,B,-A(1:end-1)./2,-A];

%Dark count in count per second
[dc, ~] = get_DC(dcname,0);

%number of pixels of the detector
n_pixl = 23;

%number of pixels in recorded image
s_pixl_x = head.ImgHdr_PixX;
s_pixl_y = head.ImgHdr_PixY;


%Spad spatial resolution distan between PM pixels in µm
SPAD_R = 23;

%scan resolution step per pixel in µm
IM_R = head.ImgHdr_PixResol;
% IM_R = 50.6826 *1e-3;

%ISM scale factor 
alpha = (1 + 523/800);

ind = im_posy<=s_pixl_y;
im_chan = im_chan(ind);
im_posx = im_posx(ind);
im_posy = im_posy(ind);


%% 
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

figure
axis equal
hold on
quiver(xc, yc ,sv_pc(1,:),sv_pc(2,:),0, 'LineWidth', 2)  
numb = 0:22;
txt = string(numb);
plot(xc, yc,'xb','MarkerSize',10, 'LineWidth', 2 ,'MarkerEdgeColor', 'green')
text(xc+0.2, yc, txt)
set(gca,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1]);
ylabel('y shift [pixel]')
xlabel('x shift [pixel]')
grid on
box on
title("SPAD pc")

xc  = -1.*mean(shift_x_ic);
yc  = -1.*mean(shift_y_ic);
% this is the shift vector for each channel
sv_ic  = -[xc; yc]./alpha;

h = figure;
ax = axes(h);
axis equal
hold on
quiver(xc, yc ,sv_ic(1,:),sv_ic(2,:),0, 'LineWidth', 2)  
numb = 0:22;
txt = string(numb);
plot(xc, yc,'xb','MarkerSize',10, 'LineWidth', 2, 'MarkerEdgeColor', 'green')
text(xc+0.2, yc, txt)

% plot(det_x, det_y,'xb','MarkerSize',10, 'LineWidth', 2, 'MarkerEdgeColor', 'red')
% text(det_x+0.2, det_y, txt)

set(ax,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1]);
ylabel('y shift [pixel]')
xlabel('x shift [pixel]')
grid on
box on
title("SPAD ic")
name = append(img_name, '_sv');
file_name = append(name,'.png');
exportgraphics(ax, file_name,'Resolution',600)

figure
axis equal
hold on
quiver(zeros(1,23), zeros(1,23) ,-sv_ic(1,:),-sv_ic(2,:),0, 'LineWidth', 2)  
set(gca,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1]);
ylabel('y shift [pixel]')
xlabel('x shift [pixel]')
grid on
box on
title("SPAD ic from zero")

M_new = zeros(n_pixl,1);
for i = 1:n_pixl
        %shift in the image
        delta_r = IM_R*sqrt( (sv_ic(1,i).*alpha)^2 + (sv_ic(2,i).*alpha)^2);
        
        de_y_i = det_x(i);
        de_x_i = det_y(i);
        
        de_y_j = det_x(11);
        de_x_j = det_y(11);
        %shift on the detector
        delta_R = SPAD_R*sqrt((de_y_i - de_y_j)^2 + (de_x_i - de_x_j)^2);

        if delta_r ~= 0
            M_new(i) = delta_R/delta_r;
        end
end

M = mean(M_new(M_new ~= 0));
disp(M/alpha)

pixel_img   = mHist(im_chan,0:22);
hex_plot(pixel_img,hot)

shift_vector = sv_ic;
save('shift_vector.m', 'shift_vector');
end

