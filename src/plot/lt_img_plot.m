function lt_img_plot(lt_img, in_img,c_map, threshold, lims, titel_name, name, sb_lenght, IM_R, pix_bin, save)
%IMG_PLOT plot 2D image with color bar and scale and save to dir
gamma = 1;

lt_img = lt_img(5:end-4,5:end-4);
in_img = in_img(5:end-4,5:end-4);

TNTvisualizer({in_img, lt_img}, struct('title',titel_name,'metadata',struct('pixelsize',IM_R,'pixelsize_unit',[char(181) 'm'])));

im = LT_image(lt_img, in_img, c_map, threshold, gamma, lims);

% %plot image 
% %reso is off here
% img_plot(im, c_map, name, titel_name, sb_lenght, IM_R, pix_bin, 0, [], 0, save, [], [])
h = figure;
ax = axes(h);

y_pix = size(lt_img,1);
x_pix = size(lt_img,2);

%position of scale bar
y1 = 0.05*y_pix;
x1 = 0.05 * x_pix;
x2 = y1 + sb_lenght/IM_R;
txt = [num2str(sb_lenght), ' µm'];
txt_offest = 0.05 * y_pix;

image(im)
colormap(c_map)
caxis(lims)
hold on;
axis off;
line([x1,x2],[y1,y1],'Color','w', 'LineWidth', 5);
text((x1+x2)/2, y1+txt_offest, txt, 'HorizontalAlignment','center', 'Color', 'w', 'FontSize', 8)
set(ax,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1], ...
    'XDir','normal', ...
    'YDir','reverse');
c = colorbar;
c.Label.String = 'Lifetime [ns]';
set(findall(h,'-property','FontSize'),'FontSize',13)
% set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)
title(titel_name);

if save
    file_name = append(name,'.png');
    exportgraphics(ax, file_name,'Resolution',600)
    spath = append(pwd,'\',name, '.tiff');
    bfsave(im, spath);
end
end


