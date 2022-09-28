function [lt_img] = get_single_lifetime(img,lin,options)
%GET_SINGLE_LIFETIME Calculates single exponential decay lifetime image via
%MLE pattern matching from single phoon data

% im_tcspc,tail_t,tail_bin,tcspc_start,max_lt,lt_cut_off
arguments
    img (:,:) double = checkerboard(256);
    lin (:,:) double = ones(1e5,1)
    %TCSPC channel of photons
    options.im_tcspc
    %edge time points of tail in ns
    options.tail_t (1,:) = ones(1,20);
    %length of tail bins in ns
    options.tail_bin = 0.05;
    %tail start time in ns
    options.tail_start_time = 3;
    %tcspc bins for tail decay
    options.tcspc_t = 1:100;
    %Maximum liftime for fit
    options.max_lt = 8;
    %Minimum number fo photons for fit
    options.lt_cut_off = 25;
end

interval    = ceil(1.01* max(max(img)));
n_pixel     = numel(img(:));
img_size    = size(img);
lt          = zeros(n_pixel, numel(options.tail_t));
n_events    = numel(lin);

h = waitbar(0,'binning');
%find all photons in one ISM pixel, by seaching in an interval of max count
%lengh + 1 
lower = 1;
upper = interval;

% Sort for faster seache
[lin, sort_index]   = sort(lin);
options.im_tcspc    = options.im_tcspc(sort_index);

for i = 1:n_pixel
   [y,x]    = ind2sub(img_size,i);
   ind      = find(lin(lower:upper) == i);
   if ~isempty(ind)
       if(img(y,x) >= options.lt_cut_off)
           [count, ~]   = histcounts(options.im_tcspc(lower + ind-1), options.tcspc_t);
           lt(i,:)      = count./sum(count);
       end
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
fprintf('lifetime pattern matching \n');

%tail fit via pattern matching
[lt_img,~]  = lt_patternMatching(lt, options.tail_t-options.tail_start_time,...
    options.tail_bin, options.max_lt);

%fit lin index to 2d array
lt_img  = reshape(lt_img, img_size);

end