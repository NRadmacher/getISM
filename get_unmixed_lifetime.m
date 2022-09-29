function [lt_amp_img] = get_unmixed_lifetime(img,lin,options)
%GET_UNMIXED_LIFETIME Calculates single exponential decay lifetime image via
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
    %name of saved image
    options.name = 'lifetime'
end

interval    = ceil(1.01* max(max(img)));
n_pixel     = numel(img(:));
img_size    = size(img);
lt_amp      = zeros(n_pixel, 4);
n_events    = numel(lin);

lt_short_name   = compose('%s lt short %0.1f ns', options.name, pattern_tau(1)); %1.37
lt_middle_name  = compose('%s lt middle %0.1f ns', options.name, pattern_tau(2)); %2.37
lt_long_name    = compose('%s lt long %0.1f ns', options.name, pattern_tau(3)); % 3.0

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
           lt_amp(i,:) = lsqnonneg(options.pattern, count.');
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

lt_amp_img = reshape(lt_amp, [img_size(1), img_size(2), 4]);

img_plot(lt_amp_img(:,:,1), c_blue, lt_short_name{1}, 'Cy2: SYT 1', 4, IM_R, reso_line_ISM, 0, 1);
img_plot(lt_amp_img(:,:,2), c_red, lt_middle_name{1}, 'OG: PSD95', 4, IM_R, reso_line_ISM, 0, 1);
img_plot(lt_amp_img(:,:,3), c_green, lt_long_name{1}, 'Alexa: GFAP', 4, IM_R, reso_line_ISM, 0, 1);

lt_rgb(:,:,1) = mat2gray(lt_amp_img(:,:,2),[2 173]);%red
lt_rgb(:,:,2) = mat2gray(lt_amp_img(:,:,3),[2 72]);%green
lt_rgb(:,:,3) = mat2gray(lt_amp_img(:,:,1),[1 73]);%blue

img_plot(lt_rgb(:,:,3), c_blue, lt_short_name{1}, 'Cy2: SYT 1', 4, IM_R, reso_line_ISM, 0, 0);
img_plot(lt_rgb(:,:,1), c_red, lt_middle_name{1}, 'OG: PSD95', 4, IM_R, reso_line_ISM, 0, 0);
img_plot(lt_rgb(:,:,2), c_green, lt_long_name{1}, 'Alexa: GFAP', 4, IM_R, reso_line_ISM, 0, 0);



end