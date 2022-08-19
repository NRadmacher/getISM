function [SOFI_img, ISM_SOFI_img] = get_SOFI(sum_img, sum_lin, im_time, im_chan, sv, head)
%GET_SOFI Calculate scanning SOFI image from single photon data and
%performs ISM pixel reasigment.

% SOFI_img      Image after pixelwise SOFI analysis
% ISM_SOFI_img  Image with addition ISM resolution enhancemnt

    %% SOFI Params

    % sum img size
    sum_size = size(sum_img);

    %number of pixels in recorded image
    s_pixl_x    = head.ImgHdr_PixX;
    s_pixl_y    = head.ImgHdr_PixY;

    % number of detector
    n_pixl = 23;

    %number of events
    n_events = numel(sum_lin);

    % SOFI shift vectors
    sofi_sv = round(sv);

    %aquisition time per scan pixel
    IM_dwell = head.ImgHdr_DwellTime;

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
    
    fprintf('Using %g ms frame lenght. Resulting in max %g frames per pixel\n', img_d*1e3, ceil(n_frames))
    fprintf('With %g frames per batch. Resulting in %g batches per pixel\n', batch_length, n_batch)
    %% SOFI ANAYSIS
    interval        = ceil(1.01* max(max(sum_img)));
    n_pixel_sum     = prod(sum_size);
    SOFI_img        = zeros(n_pixel_sum, 1);
    SOFI_ism        = zeros(n_pixel_sum, 1);
    frame_times     = linspace(0, IM_dwell, n_frames + 1);

    h = waitbar(0,'sofiing ?');
    %find all photons in one ISM pixel, by seaching in an interval of max count
    lower = 1;
    upper = interval;

    % care y x könnten vertauscht sein
    [sum_lin, sort_index] = sort(sum_lin);

    sofi_time   = im_time(sort_index);
    sofi_chan   = im_chan(sort_index);

    for i = 1:n_pixel_sum
       [y,x]            = ind2sub(sum_size,i);
       ind              = find(sum_lin(lower:upper) == i);
       if ~isempty(ind)
           
            % exact arrival time and detector channel of photons in
            % current image pixel(frame)
            f_time = sofi_time(lower + ind - 1);% + sofi_tcspc(lower + ind-1) * time_R;
            f_chan = sofi_chan(lower + ind - 1);
            
            f_time = f_time - min(f_time);
            
            f_time = combine_photon_time(f_time, 10);
            
            detector = zeros(23, 1, n_frames);
            for k = 1:n_pixl
                ind_chan = f_chan == k - 1;
                tmp      = f_time(ind_chan);
                
                [N,~]    = histcounts(tmp,frame_times);
                
                detector(k,1,:) = N;
            end
            [sof, ~] = SOFIAnalysis(detector, 2, n_frames);
            zero_lagtime = var(detector, 0, 3);
            spad_sofi = sof + zero_lagtime;
            % TO DO sum here ?
            SOFI_img(i,:) = sum(spad_sofi, 'all');
    
            sv_x = max(min(x + sofi_sv(1,:), s_pixl_x), 1);
            sv_y = max(min(y + sofi_sv(2,:), s_pixl_y), 1);
    
            lin_shift = sub2ind(sum_size, sv_y, sv_x);
    
            SOFI_ism(lin_shift) = SOFI_ism(lin_shift) + spad_sofi;
    
            % set lower edge to last found puls 1
            n_lower      = lower + ind(end);
            % set new upper edge to lower plus interval length
            upper        = min(n_lower + interval, n_events);
       else
            % set new searche boundarys. because nothing was found keep lower
            % extend upper
            n_lower      = lower;
            upper        = min(upper + interval, n_events);
       end
       % set lower
       lower        = n_lower;

       if (mod(i,n_pixel_sum/100) == 0)
            waitbar(i/n_pixel_sum,h)
       end
    end

    close(h);

    SOFI_img  = reshape(SOFI_img, sum_size);
    ISM_SOFI_img  = reshape(SOFI_ism, sum_size);
end