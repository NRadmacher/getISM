function get_wf_shift_vectors()
%GET_WF_SHIFT_VECTORS Summary of this function goes here
    
    M = 100 * 2* 75 /100;
    alpha = 23*1e-6/(M * 50*1e-9);
    det_x = -alpha.*[-2,-1,0,1,2,-3/2,-1/2,1/2,3/2,-2,-1,0,1,2,-3/2,-1/2,1/2,3/2,-2,-1,0,1,2];
    A = ones(1,5).*sqrt(3);
    B = zeros(1,5);
    det_y = -alpha.*[A,A(1:end-1)./2,B,-A(1:end-1)./2,-A];
    
    n_pixl = 23;


    h = figure;
    ax = axes(h);
    hold on
    axis equal
%     quiver(xc.*2, yc.*2 ,sv_ic(:,1),sv_ic(:,2),0, 'LineWidth', 2)  
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
    title('Wide field sv')
    file_name = append('wf_sv','.png');
    exportgraphics(ax, file_name,'Resolution',600)

    sv_wf = -1.*[det_x; det_y];

    save('SPAD_shift_vectors_wf.m', 'sv_wf');

end