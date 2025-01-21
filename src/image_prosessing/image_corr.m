function [dx,dy] = image_corr(img_a,img_b,wd)
%IMAGE_CORR Determine shift between two images using cross/image Correlation.
%returns zeros if shift cant be determined

if isempty(img_a) || isempty(img_b)
    dx = 0;
    dy = 0;
    return
end

if nargin<3
    wd = 7;
end

if any(size(img_b)~= size(img_a))
    dx = 0;
    dy = 0;
    return
end
r = zeros(wd*2+1);

for i =-wd:wd
    for j = -wd:wd
        r(i+wd+1,j+wd+1) = sum(img_a.*circshift(img_b,[i j]),'all');
    end
end

maximum = max(max(r));
[shift_y,shift_x]=find(r==maximum);

% h = figure;
% ax = axes(h);
% colormap(h,cmap_greenFireBlue)
% imagesc(ax,r);
% hold on;
% axis off;
% set(ax,'DataAspectRatio', [1,1,1], ...
%     'PlotBoxAspectRatio',[1 1 1], ...
%     'XDir','normal', ...
%     'YDir','reverse');
% file_name = append('beads',string(maximum),'.png');
% exportgraphics(ax, file_name,'Resolution',600)
% close(h)

if ~isempty(shift_y)&& ~isempty(shift_x)

    g = [-wd wd];

    [cx,cy] = meshgrid(g(1):g(2),g(1):g(2));
    [fx,fy] = meshgrid(g(1):0.1:g(2),g(1):0.1:g(2));

    cf = interp2(cx,cy,r,fx,fy,'cubic');

    maximum = max(max(cf));
    [shift_y,shift_x]=find(cf==maximum);
    if numel(shift_x)> 1
        disp("stop")
    end
    dx = fx(shift_y(1),shift_x(1));
    dy = fy(shift_y(1),shift_x(1));
else
    dx = 0;
    dy = 0;
end

