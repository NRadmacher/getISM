function rgb_img_plot(img, name, title_n, sb_lenght, IM_R, save)
%RGB_IMG_PLOT plot 2D image from RGB matrix scale and save to dir

%cut artefacts at edges
img = img(5:end-4,5:end-4,:);

h = figure;
ax = axes(h);
%actula image size
y_pix = size(img,1);
x_pix = size(img,2);

%position of scale bar
y1 = 0.05*y_pix;
x1 = 0.05 * x_pix;
x2 = y1 + sb_lenght/IM_R;
txt = [num2str(sb_lenght), ' µm'];
txt_offest = 0.05 * y_pix;

imshow(img);
hold on;
axis off;
line(ax,[x1,x2],[y1,y1],'Color','w', 'LineWidth', 5);
text(ax,(x1+x2)/2,y1+txt_offest,txt, 'HorizontalAlignment','center', 'Color', 'w', 'FontSize', 8)
set(ax,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1], ...
    'XDir','normal', ...
    'YDir','reverse');
set(findall(h,'-property','FontSize'),'FontSize',13)
title(title_n);

if save
    file_name = append(name,'.png');
    exportgraphics(ax, file_name,'Resolution',600)
%     spath = append(pwd,'\',name, '.tiff');
%     bfsave(img, spath);
end

end