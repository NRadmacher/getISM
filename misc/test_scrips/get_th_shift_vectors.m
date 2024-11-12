function sv_wf = get_th_shift_vectors()
%GET_WF_SHIFT_VECTORS Summary of this function goes here
    
    %magnificatio
    M = 100 * 2* 75 /100;
    %pixel size in [µm]
    sPixel = 0.05;
    %detector pixel size [µm]
    sPixelD = 23;
    %distance between pixels in the image [#pixels]
    alpha = sPixelD/M/sPixel;


    det_x = -alpha.*[-2,-1,0,1,2,-3/2,-1/2,1/2,3/2,-2,-1,0,1,2,-3/2,-1/2,1/2,3/2,-2,-1,0,1,2];
    A = ones(1,5).*sqrt(3);
    B = zeros(1,5);
    det_y = -alpha.*[A,A(1:end-1)./2,B,-A(1:end-1)./2,-A];
    
    n_pixl = 23;


    h = figure;
    ax = axes(h);
    hold on
    axis equal
    quiver(det_x, det_y ,-det_x./2,-det_y./2,0, 'LineWidth', 2)  
    numb = 0:n_pixl-1;
    txt = string(numb);
    plot(det_x, det_y,'xb','MarkerSize',10, 'LineWidth', 2, 'MarkerEdgeColor', 'green')
    % plot(xc.*2, yc.*2 ,'o','MarkerSize',10, 'LineWidth', 2, 'MarkerEdgeColor', 'red')
    set(ax,'DataAspectRatio', [1,1,1], ...
        'PlotBoxAspectRatio',[1 1 1]);
    ylabel('y shift [pixel]')
    xlabel('x shift [pixel]')
    XL = get(ax, 'XLim');
    xl = XL(2) - XL(1);
    text(det_x+(0.05*xl), det_y, txt)
    grid on
    box on
    %make axis perfect square
    limits = max([abs(ax.YLim), abs(ax.XLim)], [], 'all');
    ylim( [-limits, limits] );
    xlim( [-limits, limits] );
    title('thoretical sv')

    sv_wf = -0.5.*[det_x; det_y];
end



