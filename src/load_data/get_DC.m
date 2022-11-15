function [dc, bin_dc] = get_DC(fname,plt)
% Dark Count for PIIT SPAD-Array

[im_chan,im_tcspc,~,~,~, head] = read_ISM(fname);

[dc,~] = histcounts(im_chan, 23);
%count / time in sec
dc = dc./(head.ImgHdr_PixelTime * head.ImgHdr_PixNum);
dc = dc.';
% tcsps dc
bin_dc = zeros(23,head.max_bin);
for i = 1:23
   tcspc = im_tcspc(im_chan == i-1);
   [bin_dc(i,:),~] = histcounts(tcspc, 1: head.max_bin+1);
end

%count / time in sec
bin_dc = bin_dc./(head.ImgHdr_PixelTime * head.ImgHdr_PixNum);

if(plt)
    figure
    bar(0:22,dc);
    title("SPAD Darkcount")
    xlabel(sprintf('Pixel index'));
    ylabel(sprintf('counts / s'));
end

end







