function [dc, bin_dc] = get_DC(fname,plt)
% Dark Count for PIIT SPAD-Array

% n_pixl = 32;
[im_chan,im_tcspc,~,~,~, head] = read_ISM(fname);

if(sum(head.HWInpChan_Enabled,"all") > 23)
    n_pixl = 32;
    title_name = 'MPMT Darkcount';
else
    n_pixl = 23;
    title_name = 'SPAD Darkcount';
end

[dc,~] = histcounts(im_chan, n_pixl);
%count / time in sec
dc = dc./(head.ImgHdr_PixelTime * head.ImgHdr_PixNum);
dc = dc.';
mean_dc = mean(dc(dc<=max(cd,[],'all')));
% tcsps dc
bin_dc = zeros(n_pixl,head.max_bin);
for i = 1:n_pixl
   tcspc = im_tcspc(im_chan == i-1);
   [bin_dc(i,:),~] = histcounts(tcspc, 1: head.max_bin+1);
end

%count / time in sec
bin_dc = bin_dc./(head.ImgHdr_PixelTime * head.ImgHdr_PixNum);

if(plt)
    f = figure;
    ax = axes(f);
    bar(0:n_pixl-1,dc);
    yline(ax, mean_dc,'-','mean w/o hot pixel', 'Color', 'r','LineWidth',3);
    title(title_name)
    xlabel(sprintf('Pixel index'));
    ylabel(sprintf('counts / s'));
    set(ax,'YScale', 'log')
%     ylim(ax,[1 10^(ceil(log10(max(dc(:)))))])

%     get percentile

    P = prctile(dc,1:100);

    h = figure;
    ax = axes(h);
    plot(ax, 1:100, P);

    if n_pixl == 32
        pixel_img               = dc.';
        pixel_img               = log10(pixel_img);
        hex_img                 = zeros(37,1);
        hex_img(fiber_code()+1) = pixel_img;
    else
        hex_img = dc;
    end

    hex_plot(hex_img,hot)
end

end







