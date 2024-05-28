function [fwhm, indMax, t_count, im_chan, im_tcspc] = get_IRF(fname, plt)
%% Load Data

if (nargin<1)
   %TO DO find way to save and fast access! 
end
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220204\irf_2ph_014.ptu';
% fname = 'D:\PHD\Data\2022\220816\irf_ex470nm_004.ptu';
% dcname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\210909\IRF_DC_001.ptu';
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230130_pmt_irf\irf_pmt_011.ptu';

[im_chan,im_tcspc,~,head] = read_FCS(fname);

% im_tcspc = remove_MHH_offset(im_tcspc,im_chan, head.max_bin, 500);

bin_factor = 1;

if(nnz(head.TTResult_InputRate) > 23)
    title_name = 'MPMT irf';
else
    title_name = 'SPAD-Array IRF';
end

%% find fwhm and max
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
    legende_txt = append('IRF', newline, 'FWHM: ', string(fwhm*time_R*1e12), ' ps');
    legend(legende_txt, 'Location', 'northeast')
    grid(ax, 'on')
    xlabel(sprintf('time [ns]'));
    ylabel(sprintf('count'));
    xlim([3 11])
    title(title_name)
    set(ax, 'FontSize', 13, 'FontWeight', 'bold', 'YScale', 'log', ...
        'Box', 'on', 'LineWidth', 1)

    figure
    hold on
    for i = 1:32
        [c_count, ~] = histcounts(im_tcspc(im_chan == i-1), 1:bin_factor:max_bin+1);
        chr = int2str(i-1);
        plot(time,c_count, 'DisplayName',chr)
    end
    title(img_name);

    file_name = append(img_name,'.png');
    exportgraphics(ax, file_name,'Resolution',600)
end

%save irf 
% count = count.';
% save('irf.m', 'count');

end