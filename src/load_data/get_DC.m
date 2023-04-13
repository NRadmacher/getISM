function [dc, bin_dc] = get_DC(fname,plt)
% Dark Count for PIIT SPAD-Array

% n_pixl = 32;
[im_chan,im_tcspc,~,~,~, head] = read_ISM(fname);

if(sum(head.HWInpChan_Enabled,"all") > 23)
    n_pixl = 32;
    title_name = 'MPMT Darkcount';
else
    n_pixl = 23;
    title_name = 'SPAD-Array Darkcount';
end

[dc,~] = histcounts(im_chan, n_pixl);
%count / time in sec
[~,ind] = max(dc);
% dc(ind) = 0;
dc = dc./(head.ImgHdr_PixelTime * head.ImgHdr_PixNum);
mean_dc = median(dc);
% tcsps dc
bin_dc = zeros(n_pixl,head.max_bin);
for i = 1:n_pixl
   tcspc = im_tcspc(im_chan == i-1);
   [bin_dc(i,:),~] = histcounts(tcspc, 1: head.max_bin+1);
end

%count / time in sec
bin_dc = bin_dc./(head.ImgHdr_PixelTime * head.ImgHdr_PixNum);

if(plt)

    tmp_name = strsplit(fname, '\');
    date = tmp_name{end-1};
    img_name = tmp_name{end};
    img_name = strsplit(img_name, '.');
    img_name = img_name{end-1};
    img_name = append(date,' ',img_name);
    img_name = strrep(img_name,'_',' ');

    f = figure;
    ax = axes(f);
    bar(0:n_pixl-1,dc,'FaceColor',[0.00,1.00,0.00],'EdgeColor', [0.39,0.83,0.07],'FaceAlpha',1);
    line_txt = append('median: ', string(mean_dc));
    yline(ax, mean_dc,'-',line_txt, 'Color', 'r','LineWidth',3, 'Alpha',1,'FontSize',13,'FontWeight','bold');
    title(title_name)
    xlabel(sprintf('Pixel index'));
    ylabel(sprintf('counts / s'));
    set(ax,'YScale', 'log', 'FontSize',13,'FontWeight','bold', 'LineWidth',1)
    grid(ax,"on")
    ylim(ax,[10^(floor(log10(min(dc(:))))) 10^(ceil(log10(max(dc(:)))))])
    ylim(ax,[1 10^(ceil(log10(max(dc(:)))))])

    file_name = append(img_name,'_bar.png');
    exportgraphics(ax, file_name,'Resolution',600)

%     get percentile

    P = prctile(dc,1:100);

    h = figure;
    ax = axes(h);
    plot(ax, 1:100, P,'LineWidth',2);
    xlabel(sprintf('Percentage of Pixels'));
    ylabel(sprintf('counts / s'));
    title(title_name)
    grid(ax,"on")
    set(ax,'YScale', 'log', 'FontSize',13,'FontWeight','bold', 'LineWidth',1)

    file_name = append(img_name,'_per.png');
    exportgraphics(ax, file_name,'Resolution',600)
    
    if n_pixl == 32
        pixel_img               = dc.';
        pixel_img               = (pixel_img);
        hex_img                 = zeros(37,1);
        hex_img(fiber_code()+1) = pixel_img;

        pmt_plot(im_chan, n_pixl, hot)
    else
        hex_img = dc;
    end

    hex_plot(hex_img, hot)
end

end







