function [lt_amp_img] = get_bi_lifetime(img,lin,options)
%GET_BI_LIFETIME Calculates single exponential decay lifetime image via
%MLE pattern matching from single phoon data

arguments
    img (:,:) double = checkerboard(256);
    lin (:,1) double = ones(1e5,1)
    %TCSPC channel of photons
    options.im_tcspc (:,1) double = ones(1e5,1)
    %edge time points of tail in ns
    options.tail_t (1,:) = ones(1,20);
    %length of tail bins in ns
    options.tail_bin = 0.05;
    %tail start time in ns
    options.tail_start_time = 3;
    %tcspc bins for tail decay
    options.tcspc_t = 1:100;
    %Pre calculated decay patterns for unmixing
    options.pattern
    %filetimes corresponding to patterns
    options.pattern_tau
    %Struckture and Fluorophore name
    options.name (2,2) cell
    %file name
    options.fname
    %plot of analysis
    options.ana_plt = false
end

interval    = ceil(1.01* max(max(img)));
n_pixel     = numel(img(:));
img_size    = size(img);
lt_amp      = zeros(n_pixel, 3);
n_events    = numel(lin);

if(options.ana_plt)
    
    % tail time that starts at 0
    xx = options.tail_t-options.tail_start_time;

    % tcspc histogram
    [count, ~]   = histcounts(options.im_tcspc, options.tcspc_t);
    
    % find optimal doulbne exp decy
    [double_amp,tau] = double_exp(xx.',count.');
    disp('over all tau and amp from fit')
    disp(tau)
    disp(double_amp./sum(double_amp))
    
    % find amplidtude for patterns
    pattern_amp = lsqnonneg(options.pattern, count.');
    disp('over all tau and amp from patterns')
    disp(options.pattern_tau)
    disp(pattern_amp.'./sum(pattern_amp))

    % plot decay and fits
    h = figure;
    ax = axes(h);
    
    p1 = plot(ax,xx, count, 'o', 'MarkerSize', 10, 'Color', '#EDB120');
    hold on

    amp1 = double_amp(1);
    amp2 = double_amp(2);
    offset = double_amp(3);
    y2 =  amp1*exp(-xx/tau(1)) + amp2*exp(-xx/tau(2)) + offset;
    p2 = plot(ax, xx, y2, '-', 'LineWidth', 2, 'Color', 'blue');
    
    amp1 = pattern_amp(1);
    amp2 = pattern_amp(2);
    offset = pattern_amp(3);
    y3 =  options.pattern(:,1).*amp1 + options.pattern(:,2).*amp2 + options.pattern(:,3).*offset;
    p3 = plot(ax, xx, y3, '-', 'LineWidth', 2, 'Color', '#A2142F');

    legend('tail decay','fitt', 'pattern')
    set(ax, 'YScale', 'log')
    xlim('padded')
    ylim('padded')
    xlabel('time [ns]')
    ylabel('count')
    title('TCSPC decay with fit');
    set(findall(h,'-property','FontSize'),'FontSize',15)
end

h = waitbar(0,'pixel wise pattern matching');
%find all photons in one ISM pixel, by seaching in an interval of max count
%lengh + 1 
lower = 1;
upper = interval;

% Sort for faster seache
[lin, sort_index]   = sort(lin);
tcspc    = options.im_tcspc(sort_index);

for i = 1:n_pixel
   ind      = find(lin(lower:upper) == i);
   if ~isempty(ind)
       [count, ~]   = histcounts(tcspc(lower + ind-1), options.tcspc_t);
       lt_amp(i,:) = lsqnonneg(options.pattern, count.');
       % set new upper edge to lower plus interval length
       upper    = min(lower + ind(end) + interval, n_events);
       % set lower edge to last found puls 1
       lower    = lower + ind(end);
   else
       % set new searche boundarys. because nothing was found keep
       % lower extend upper
       upper    = min(upper + interval, n_events);
   end
   if (mod(i,n_pixel/100) < 1)
        waitbar(i/n_pixel,h)
   end
end
close(h);

lt_amp_img = reshape(lt_amp, [img_size(1), img_size(2), 3]);
end
