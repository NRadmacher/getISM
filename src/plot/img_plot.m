function img_plot(img, c_map, name, title_n, sb_lenght, IM_R, pix_bin, reso_line, reso, save)
%IMG_PLOT plot 2D image with color bar and scale and save to dir

% roi = [50 1; 124 75];
% img = img( roi(1,2):roi(2,2), roi(1,1) : roi(2,1));

%cut artefacts at edges
img = img(3:end-3,3:end-3);

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
set(findall(h,'-property','FontSize'),'FontSize',13)
title(title_n);

if reso
    
    xb1 = reso_line(1,1);
    xb2 = reso_line(1,2);
    yb1 = reso_line(2,1);
    yb2 = reso_line(2,2);
    
    %horizontal limits for intensity line
    xb1 = xb1 - roi(1,1);
    xb2 = xb2 - roi(1,1);
    dxb = xb2 - xb1;
    
    %vertical limits of intensity line
    yb1 = yb1 - roi(1,2);
    yb2 = yb2 - roi(1,2);
    dyb = yb2 - yb1;
    
%     Set long and short "differnec" for sampeling along intensity line
    if(dxb > dyb)
        lb1 = xb1;
        lb2 = xb2;
        dlb = dxb;
        sb1 = yb1;
        sb2 = yb2;
        dsb = dxb;
    else
        lb1 = yb1;
        lb2 = yb2;
        dlb = dyb;
        sb1 = xb1;
        sb2 = xb2;
        dsb = dxb;
    end
    %finde coordinats along intensity line
    if(dsb ==0) %take care of vertical and horzontal lines
        llb     = lb1:1:lb2;
        ssb     = zeros(1,dlb + 1) + sb1;
        lb      = IM_R;
        r_bead  = (0:dlb) * lb;
    else %itensity line is  hypotenuse llb and ssb are the legs
        llb = lb1:1:lb2;
        ssb = ceil( dsb/dlb * (0:dlb) + sb1);
        %lenght of eatch segmant along intesity line
        lb = sqrt(dlb^2 + dsb^2) /dlb * IM_R;
        r_bead  = (0:dlb) * lb;
    end
    
    lin = sub2ind(size(img), llb, ssb);
    bead_sum = img(lin);
    bead_sum = bead_sum;%/max(bead_sum);

    fitt = fit(r_bead.', bead_sum.', 'gauss1');

    sig = fitt.c1/2;
    ci = confint(fitt,0.95)/2;
    sig_err = abs(sig - max(ci(:,3)));

%     fwhm = 2*sqrt(2 * log(2)) * sig;
%     fwhm_err = 2*sqrt(2 * log(2)) * sig_err;
    
    fwhm = 2 * sig;
    fwhm_err = 2 * sig_err;

    %get significan decimal of error
    sig_deci = log10(fwhm_err);
    sig_deci = abs(round(sig_deci)) + 1;

    fwhm_err = round(fwhm_err, sig_deci);
    fwhm = round(fwhm, sig_deci);

    fit_stg = sprintf(' %g \x00B1 %g nm', fwhm*1000, fwhm_err*1000);


    r = figure; 
    plot(fitt, r_bead, bead_sum,'bx')
    set(findall(r,'-property','FontSize'),'FontSize',17)
    set(findall(r,'-property', 'MarkerSize'), 'MarkerSize', 12)
    set(findall(r,'-property','LineWidth'),'LineWidth',1.5)
%     ylim([0 1.1])
    xlabel(sprintf('distance / %sm',char(181)));
    ylabel('intensity');
    legend({' data', [' gaussian fit' newline ' resolution:' newline fit_stg]}, 'Location', 'northwest')
    grid on
    line(ax,[xb1,xb2],[yb1,yb2],'Color','r', 'LineWidth', 2);
end

if save
    file_name = append(name,'.png');
    exportgraphics(ax, file_name,'Resolution',600)
    spath = append(pwd,'\',name, '.tiff');

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

