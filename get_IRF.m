function [fwhm, indMax, t_count, im_chan, im_tcspc] = get_IRF(fname, plt)
%% Load Data

if (nargin<1)
   %TO DO find way to save and fast access! 
end
% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220204\irf_2ph_014.ptu';
% % fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220106\irf_007.ptu';
% dcname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\210909\IRF_DC_001.ptu';

% [im_chan,im_tcspc,~,~,~,head] = read_ISM(fname);
[im_chan,im_tcspc,~,head] = read_FCS(fname);

bin_factor = 1;

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
    
    tmp_name = strsplit(fname, '\');
    date = tmp_name{end-1};
    img_name = tmp_name{end};
    img_name = strsplit(img_name, '.');
    img_name = img_name{end-1};
    img_name = append(date,' ',img_name);
    img_name = strrep(img_name,'_',' ');
    
    time_R = mean(head.MeasDesc_Resolution);
    figure
    hold on
    plot(count, '.')
    set(gca, 'YScale', 'log')

    time = 1:bin_factor:max_bin;
    time = time.* (time_R *1e9);
    figure
    hold on
    for i = 1:23
        [c_count, ~] = histcounts(im_tcspc(im_chan == i-1), 1:bin_factor:max_bin+1);
    %     c_count = c_count/max(c_count);
        chr = int2str(i-1);
        plot(time,c_count, 'DisplayName',chr)
    end
    % xlim([0 3])
    title(img_name);
    legend
end

%save irf 
% count = count.';
% save('irf.m', 'count');

end