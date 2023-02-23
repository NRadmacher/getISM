function [fwhm, indMax, t_count, im_chan, im_tcspc] = get_IRF(fname, plt)
%% Load Data

if (nargin<1)
   %TO DO find way to save and fast access! 
end
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220204\irf_2ph_014.ptu';
% fname = 'D:\PHD\Data\2022\220816\irf_ex470nm_004.ptu';
% dcname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\210909\IRF_DC_001.ptu';
fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230130\irf_pmt_011.ptu';

[im_chan,im_tcspc,~,head] = read_FCS(fname);

im_tcspc = remove_MHH_offset(im_tcspc,im_chan, head.max_bin, 500);

bin_factor = 1;

if(sum(head.HWInpChan_Enabled,"all") > 23)
    title_name = 'MPMT';
else
    title_name = 'SPAD-Array';
end

%% find fwhm and max
% max_bin = round(1/head.TTResult_SyncRate /head.HW_BaseResolution);
max_bin = head.max_bin;

[count, ~] = histcounts(im_tcspc, 1:bin_factor:max_bin+1);
% Find the half max value.
[fullMax, indMax] = max(count);
halfMax = (min(count) + fullMax) / 2;
% Find where the data first drops below half the max.
index1 = find(count >= halfMax, 1, 'first');
% Find where the data last rises above half the max.
index2 = find(count >= halfMax, 1, 'last');
fwhm = index2-index1 + 1; % FWHM in indexes.
t_count = numel(im_tcspc);

fprintf('FWHM: %d, Max@ %d\n',fwhm, indMax );
%% plot sum irf and for eatch pixel 
if(plt)
    time_R = mean(head.MeasDesc_Resolution);
    time = 1:bin_factor:max_bin;
    time = time.* (time_R *1e9);
    
    tmp_name = strsplit(fname, '\');
    date = tmp_name{end-1};
    img_name = tmp_name{end};
    img_name = strsplit(img_name, '.');
    img_name = img_name{end-1};
    img_name = append(date,' ',img_name);
    img_name = strrep(img_name,'_',' ');
    
    h = figure;
    ax = axes(h);
    hold on
    plot(time,count, '.','LineStyle','-','Marker','none',LineWidth=1)
    set(ax, 'FontSize', 13, 'FontWeight', 'bold', 'YScale', 'log')
    legende_txt = append('IRF', newline, 'FWHM: ', string(fwhm*time_R*1e12), ' ps');
    legend(legende_txt, 'Location', 'northeast')
    grid('on')
    xlabel(sprintf('time [ns]'));
    ylabel(sprintf('count'));
    xlim([0 10])
    title(title_name)

    figure
    hold on
    for i = 1:32
%         if(i == 11 || i == 15)
        [c_count, ~] = histcounts(im_tcspc(im_chan == i-1), 1:bin_factor:max_bin+1);
    %     c_count = c_count/max(c_count);
        chr = int2str(i-1);
        plot(time,c_count, 'DisplayName',chr)
%         end
    end
    % xlim([0 3])
    title(img_name);
end

%save irf 
% count = count.';
% save('irf.m', 'count');

end