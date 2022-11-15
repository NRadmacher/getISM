function ISM_SOFI_img = get_ISMSOFI(ism_img, ism_lin, im_time, head)
%GET_SOFI Calculate scanning SOFI image from single photon data and
%performs ISM pixel reasigment.

% SOFI_img      Image after pixelwise SOFI analysis
% ISM_SOFI_img  Image with addition ISM resolution enhancemnt

    %% SOFI Params

    % sum img size
    ism_size = size(ism_img);

    %number of events
    n_events = numel(ism_lin);

    %aquisition time per scan pixel
    IM_dwell = head.ImgHdr_PixelTime;

    pix_dwell = head.ImgHdr_DwellTime;

    %duration of one frames (image) in seconds for SOFI JE 100 mu sec
    img_d = 1e-4;
    
    %total number for frames per pixel
    n_frames = round(IM_dwell /img_d);
    
    %number of frames per batch
    batch_length = 500;
    
    % fix ?
    if n_frames < batch_length
        fprintf('more frames per batch than total frames per pixel. Setting batch size to n_frames\n')
        batch_length = n_frames;
    end
    
    %number of bathes with batch_length images
    n_batch = ceil(n_frames / batch_length);
    
    %fix ?
    if(n_batch * batch_length > n_frames)
        n_batch = n_batch - 1;
    end
    
    fprintf('Using %g ms frame lenght. Resulting in max %g frames per pixel. With %g s total dwell time\n', img_d*1e3, ceil(n_frames), IM_dwell)
    fprintf('With %g frames per batch. Resulting in %g batches per pixel\n', batch_length, n_batch)
    %% SOFI ANAYSIS
    
    interval        = ceil(1.01* max(max(ism_img)));
    n_pixel_ism     = prod(ism_size);
    SOFI_ism        = zeros(n_pixel_ism, n_frames);
    SOFI_time       = zeros(n_pixel_ism, 1);
    frame_times     = linspace(0, IM_dwell, n_frames + 1);

    %find all photons in one ISM pixel, by seaching in an interval of max count
    lower = 1;
    upper = interval;

    % care y x könnten vertauscht sein
    % problem with ISM than SOFI or in other words doing SOIF on ISM ist
    % the timeming problem of photons. ISM pixel gets photons fom many scan
    % pixel
    [ism_lin, sort_index] = sort(ism_lin);

    sofi_time   = im_time(sort_index);

    h = waitbar(0,'sofiing ?');
    % going over all pixels
    for i = 1:n_pixel_ism
       [y,x]            = ind2sub(ism_size,i);
       ind              = find(ism_lin(lower:upper) == i);
       if((x==127) && (y ==57))
            disp('stop')
       end
       if ~isempty(ind)
           
            % exact arrival time and detector channel of photons in
            % current image pixel(frame)
            f_time = sofi_time(lower + ind - 1);% + sofi_tcspc(lower + ind-1) * time_R;
            
            f_time = f_time - min(f_time);
            % remove time gaps due to multiple scans
            f_time = combine_photon_time(f_time, pix_dwell);

            [count, ~]   = histcounts(f_time, 'BinWidth',img_d);
           % SOFI_ism(i,:) = count;
            SOFI_time(i,1) = var(count);
           % set new upper edge to lower plus interval length
           upper    = min(lower + ind(end) + interval, n_events);
           % set lower edge to last found puls 1
           lower    = lower + ind(end);
       else
           % set new searche boundarys. because nothing was found keep
           % lower extend upper
           upper    = min(upper + interval, n_events);
       end
       if (mod(i,n_pixel_ism/100) < 1)
            waitbar(i/n_pixel_ism,h)
       end
    end

    close(h);

    SOFI_ism = reshape(SOFI_ism, [ism_size(1), ism_size(2), n_frames]);
    % get sum of cumulants for non zeros lagtime
    [sof, ~] = SOFIAnalysis(SOFI_ism, 2, n_frames);
    % second order cumulant with zero lagtime is just varinace
    zero_lagtime = var(SOFI_ism, 0, 3);
    ISM_SOFI_img = sof + zero_lagtime;
    ISM_SOFI_img = reshape(SOFI_time, [ism_size(1), ism_size(2)]);
end


