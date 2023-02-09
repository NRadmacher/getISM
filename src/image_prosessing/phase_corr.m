function [dx,dy] = phase_corr(img_a,img_b)
%PHASE_CORR Determine shift between two images using Phase Correlation.
%returns zeros if shift cant be determined

%scale image to the next power of two
p = max(nextpow2(size(img_a)));
extra = 2^p;
extra_h = 2^(p-1);
wd = 7;

Gb = fft2(img_b,extra,extra);
Ga = fft2(img_a,extra,extra);

R = Ga.*conj(Gb)./(abs(Ga.*conj(Gb)));
r = ifftshift(ifft2(R));

maximum = max(max(r));
[shift_y,shift_x]=find(r==maximum);

if ~isempty(shift_y)&& ~isempty(shift_x)
    %get shift relative to the center
    %there seems to be an gloabl ofset of -1 -1??
    %expression in the bracke is the shift beteween imgs
    %now invert shift

    r_sub = r(extra_h-wd : extra_h+wd, extra_h-wd : extra_h+wd);

    g = [-wd wd];

    [cx,cy] = meshgrid(g(1):g(2),g(1):g(2));
    [fx,fy] = meshgrid(g(1):0.1:g(2),g(1):0.1:g(2));

    cf = interp2(cx,cy,r_sub,fx,fy,'cubic');

    maximum = max(max(cf));
    [shift_y,shift_x]=find(cf==maximum);
    dx = fx(shift_y,shift_x);
    dy = fy(shift_y,shift_x);
else
    dx = 0;
    dy = 0;
end

