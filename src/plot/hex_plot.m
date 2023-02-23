function hex_plot(data,map)
%HEX_PLOT plots scalar data in hexagonal patches

if numel(data) == 23
    name = 'SPAD';

    %   generate centers of hexagonal lattic a = 1/chenter of the patches
    det_x = -1.*[-2,-1,0,1,2,-3/2,-1/2,1/2,3/2,-2,-1,0,1,2,-3/2,-1/2,1/2,3/2,-2,-1,0,1,2];
    det_x = repmat(det_x, 6,1);
    A = ones(1,5).*sqrt(3);
    B = zeros(1,5);
    det_y = 1.*[A,A(1:end-1)./2,B,-A(1:end-1)./2,-A];
    det_y = repmat(det_y, 6,1);
    
    numb = 0:22;

else

    name = 'MPMT';

    %   generate centers of hexagonal lattic a = 1/chenter of the patches
    det_x = 1.*[-3/2:3/2, -4/2:4/2, -5/2:5/2, -6/2:6/2, -5/2:5/2, -4/2:4/2, -3/2:3/2];
    det_x = repmat(det_x, 6,1);
    A = sqrt(3)/2;
    det_y = 1.*[3*A.*ones(1,4), 2*A.*ones(1,5), 1*A.*ones(1,6), 0*A.*ones(1,7), -1*A.*ones(1,6), -2*A.*ones(1,5), -3*A.*ones(1,4)];
    det_y = repmat(det_y, 6,1);
    
    numb = 0:36;

end

% generate shifts for the veticies of hexagonal patches
a = 1/2;
b = a * 2/sqrt(3);

y_shift = b * [1/2, 1, 1/2, -1/2, -1, -1/2]; 
x_shift = a*[1, 0, -1, -1, 0, 1];

X = bsxfun(@plus, det_x, x_shift.');
Y = bsxfun(@plus, det_y, y_shift.');
C = data;

txt = string(numb);

figure
colormap(map)
patch(X,Y,C)
text(det_x(1,:), det_y(1,:), txt)
set(gca,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1])
colorbar
title(name)
grid("off")
axis("off")
end

