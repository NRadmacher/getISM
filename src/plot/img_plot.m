function img_plot(img, c_map, name, title_n, sb_lenght, IM_R, pix_bin, roi, reso_line, reso, save, f, f_ax)
%IMG_PLOT plot 2D image with color bar and scale and save to dir

if roi == 0
    roi = [1 1; size(img,1) size(img,2)];
end
% ymin xin ;ymax xmax
img = img( roi(1,1):roi(2,1), roi(1,2) : roi(2,2));

%cut artefacts at edges
% img = img(3:end-3,3:end-3);

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

colormap(h,c_map)
imagesc(ax,img);
hold on;
axis off;
line(ax,[x1,x2],[y1,y1],'Color','w', 'LineWidth', 5);
text(ax,(x1+x2)/2,y1+txt_offest,txt, 'HorizontalAlignment','center', 'Color', 'w', 'FontSize', 8)
set(ax,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1], ...
    'XDir','normal', ...
    'YDir','reverse');
c = colorbar;
if max(img,[],'all') == 1
    c.Label.String = 'rel. intensity';
else
    c.Label.String = 'Count';
end
colorbar('off')
set(findall(h,'-property','FontSize'),'FontSize',13)
% title(append(title_n, ' image'),'FontSize',15);

if reso
    
    xb1 = reso_line(1,1);
    xb2 = reso_line(1,2);
    yb1 = reso_line(2,1);
    yb2 = reso_line(2,2);
    
    %horizontal limits for intensity line
    xb1 = xb1 - +1;
    xb2 = xb2 - +1;
    dxb = xb2 - xb1;
    
    %vertical limits of intensity line
    yb1 = yb1 - +1;
    yb2 = yb2 - +1;
    dyb = yb2 - yb1;
    
%     Set long and short "differnec" for sampeling along intensity line
%     if(dxb > dyb)
%         lb1 = xb1;
%         lb2 = xb2;
%         dlb = dxb;
%         sb1 = yb1;
%         sb2 = yb2;
%         dsb = dyb;
%     else
%         lb1 = yb1;
%         lb2 = yb2;
%         dlb = dyb;
%         sb1 = xb1;
%         sb2 = xb2;
%         dsb = dxb;
%     end
    %finde coordinats along intensity line
    if(dyb == 0 || dxb == 0) %take care of vertical and horzontal lines
        if(dyb == 0)
            xb     = xb1:1:xb2;
            yb     = zeros(1,dxb + 1) + yb1;
            lb      = IM_R/pix_bin;
            r_bead  = (0:dxb) * lb;
        else
            yb     = yb1:1:yb2;
            xb     = zeros(1,dyb + 1) + xb1;
            lb      = IM_R/pix_bin;
            r_bead  = (0:dyb) * lb;
        end
    else %itensity line is  hypotenuse llb and ssb are the legs
        llb = lb1:1:lb2;
        ssb = ceil( dsb/dlb * (0:dlb) + sb1);
        %lenght of eatch segmant along intesity line
        lb = sqrt(dlb^2 + dsb^2) /dlb * IM_R/pix_bin;
        r_bead  = (0:dlb) * lb;
    end
    
    lin = sub2ind(size(img), yb, xb);
    bead_sum = img(lin);
    bead_sum = bead_sum./max(bead_sum);
    fitt = fit(r_bead.', bead_sum.', 'gauss1');

    sig = fitt.c1/sqrt(2);
    ci = confint(fitt,0.682);
    sig_err = abs( sig - max(ci(:,3))/sqrt(2) );
    
    two_sig = 2 * sig;
    fwhm_err = 2 * sig_err;

    %get significan decimal of error
    sig_deci = log10(fwhm_err);
    sig_deci = abs(round(sig_deci)) + 1;

    fwhm_err = round(fwhm_err, sig_deci);
    two_sig = round(two_sig, sig_deci);

    fit_stg = sprintf('resolutio: %g \x00B1 %g nm', two_sig*1000, fwhm_err*1000);
    fprintf(append(name, ' ', fit_stg, '\n'));

%     r = figure;
%     r_ax = axes('Parent', r);
    r = f;
    axes(f_ax)
    hold on
    p = plot(r_bead, bead_sum, 'x', 'DisplayName', title_n);
    xx = linspace(min(r_bead), max(r_bead), 1000);
    yy = feval(fitt, xx);
    plot(xx, yy,'Color', p.Color)
%     plot(xx, yy,'Color', p.Color, 'DisplayName', fit_stg)
    set(findall(r,'-property','FontSize'),'FontSize',17)
    set(findall(r,'-property', 'MarkerSize'), 'MarkerSize', 12)
    set(findall(r,'-property','LineWidth'),'LineWidth',1.5)
%     ylim([0 1.1])
    xlabel(sprintf('distance [%sm]',char(181)));
    ylabel('rel. intensity');
%     legappend({' data', [' gaussian fit' newline ' resolution:' newline fit_stg]}, 'Location', 'northwest')
%     legend
    grid on
    line(ax,[xb1,xb2],[yb1,yb2],'Color','r', 'LineWidth', 2);
end

if save
    file_name = append(name,'.png');
    exportgraphics(ax, file_name,'Resolution',600)
   
    if reso
%         file_name = append(name,'_reso.png');
%         exportgraphics(f_ax, file_name,'Resolution',600)
    end

    spath = append(pwd,'\',name, '.tif');

    if exist(spath, 'file')==2
        delete(spath);
    end

    metadata = createMinimalOMEXMLMetadata(img);
    pix_size = IM_R/pix_bin;
    pixelSize = ome.units.quantity.Length(java.lang.Double(pix_size), ome.units.UNITS.MICROMETER);
    metadata.setPixelsPhysicalSizeX(pixelSize, 0);
    metadata.setPixelsPhysicalSizeY(pixelSize, 0);

    bfsave(img, spath,'metadata', metadata);
end

end

