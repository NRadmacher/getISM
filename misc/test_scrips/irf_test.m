clear
close all
clc


names = dir('C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220207\irf_2ph_gold_bead_*.ptu');

n_points_x = 3;
n_points_y = 3;
n_points = n_points_x * n_points_y;
n_points = 5;
start_offset = 0;
pos         = zeros(n_points,2);
irf_max     = zeros(n_points_x,n_points_y);
irf_count   = zeros(n_points_x,n_points_y);
for i = 1:n_points
    
    j = i + start_offset;
    fname = append(names(j).folder,'\', names(j).name);
    head  = PTU_Read_Head(fname);
    
    pos(i,1) = head.ImgHdr_X0;
    pos(i,2) = head.ImgHdr_Y0;
    
    [~, irf_max(i),irf_count(i), im_chan,~] = get_IRF(fname, 1);
    [pixel_img, ~]   = histcounts(im_chan,0:23);
    hex_plot(pixel_img.',hot)
end
%% plot
figure
bar3(irf_max)
xlabel('x')
xticklabels({'-75','-50','-25','0','25','50','75'})
ylabel('y')
yticklabels({'-75','-50','-25','0','25','50','75'})


figure
plot(irf_count(:),irf_max(:), 'x')



