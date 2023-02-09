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
% r = xcorr2(img_a,img_b);

maximum = max(max(r));
[shift_y,shift_x]=find(r==maximum);

if ~isempty(shift_y)&& ~isempty(shift_x)

%     r_sub = r(half_y-wd : half_y+wd, half_x-wd : half_x+wd);

    g = [-wd wd];

    [cx,cy] = meshgrid(g(1):g(2),g(1):g(2));
    [fx,fy] = meshgrid(g(1):0.1:g(2),g(1):0.1:g(2));

    cf = interp2(cx,cy,r,fx,fy,'cubic');

    maximum = max(max(cf));
    [shift_y,shift_x]=find(cf==maximum);
    dx = fx(shift_y,shift_x);
    dy = fy(shift_y,shift_x);
%     dx = shift_x - size(img_a,2);
%     dy = shift_y - size(img_a,1);
else
    dx = 0;
    dy = 0;
end

