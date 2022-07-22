%%
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F)
clear
close all
clc


lifetime    = 0;
deconv      = 0;
sofi        = 1;
%% Loadind data
% load('scimaps.mat');
%rgb values for color map black,blue,cyan,green,yellow,orange?,red,magenta
sp1     = 1:255/7:256;
sp2     = [[0 0 0]; [0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]; [1 0 1]];
br1     = 1:255/5:256;
br2     = [[0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]];
grl     = 1:255/3:256;
gr2     = [[0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]];
gl      = 1:255/1:256;
g2      = [[0 0 0]; [0 1 0]];

lambda  = 1:256;

spectrum    = interp1(sp1, sp2, lambda);
c_greenred  = interp1(grl, gr2, lambda);
c_bluered   = interp1(br1, br2, lambda);
c_green     = interp1(gl, g2, lambda);
c_map       = cmap_isoluminant75;

fname   = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220722\qdots_em585nm_007.ptu';
dcname  = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220106\cd_001.ptu';
irfname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220708\irf_ex470nm_005.ptu';

%Image title and name from file name
tmp_name        = strsplit(fname, '\');
date            = tmp_name{end-1};
img_name        = tmp_name{end};
img_name        = strsplit(img_name, '.');
img_name        = img_name{end-1};
img_name        = append(date,' ',img_name);
ISM_name        = append(img_name, ' ISM');
ISM_docn_name   = append(img_name, ' ISM Decon');
colf_name       = append(img_name, ' no ISM');
LT_name         = append(img_name, ' lt');

%read ISM data from .ptu file
[im_chan,im_tcspc,im_posy,im_posx,im_time, head] = read_ISM(fname);

ind             = im_posy<=head.ImgHdr_PixY;
im_chan         = im_chan(ind);
im_posx         = im_posx(ind);
im_posy         = im_posy(ind);
im_posx         = im_posx * head.ImgHdr_PixX;
im_pix          = im_chan;
im_time         = im_time./head.TTResult_SyncRate; % photon arrival in seconds

%remove time delay due to unsyncronised Multi Harps
% im_tcspc = remove_MHH_offset(im_tcspc,im_chan, head.max_bin, 100);

%Dark count in count per second
[dc, bin_dc] = get_DC(dcname,0);

%get IRF to set tail for tailfit
% [fwhm, indMax, ~, ~,irf_tcspc] = get_IRF(irfname, 1);
%% Parameters and magic numbers(please fix)
%number of pixels of the detector
n_pixl = 23;

%number of pixels in recorded image
s_pixl_x    = head.ImgHdr_PixX;
s_pixl_y    = head.ImgHdr_PixY;
s_pixl      = s_pixl_x * s_pixl_y;
s_size      = [s_pixl_x s_pixl_y];

%statrting position of the scan in the sample
s_start_x   = head.ImgHdr_X0;
s_start_y   = head.ImgHdr_Y0;

%Spad spatial resolution distan between PM pixels in µm
SPAD_R = 23;

%scan resolution step per pixel in µm
IM_R = head.ImgHdr_PixResol;

%aquisition time per scan pixel
IM_dwell = head.ImgHdr_DwellTime;

%ISM Binning UGO[1.02] PQ[1.04] UGO220302[1.01] 20nm[1.0075]
ISM_binning = 1.0075;

%temporal resolution of TCSPC in seconds
time_R = mean(head.MeasDesc_Resolution);

%Time between laser exciation pulses in seconds
rep_time = 1/head.TTResult_SyncRate;

%max number of tcspc bins %%CARE
max_bin = head.max_bin;

%Magnification at detector
M = 200;

%interval lenght to find tscpc times corresponding to the same image Pixel
%!!FIX NEEDED!!
interval = 8000;

% TCSPC binning factor
bin_factor = 10;

%minimum numbers of photons to calculate lifetime
lt_cut_off = 350;

% tcspc tail_start in tcspc bin number !!FIX NEEDED!!
% tail_start = indMax + fwhm;
tail_start = 1550;
% tail_start = 270;

% tcspc tail_end in tcspc bin number
tail_end = max_bin - 100;

% TCSPC bin length in ns
tcspc_bin_l = time_R *1e9;

% combin binning for lifetime fitting
% Lifetime bin lenght in ns
lt_bin_l = tcspc_bin_l * bin_factor;

% Lifetime tail 
lt_start    = floor(tail_start / bin_factor);
lt_end      = floor(tail_end / bin_factor); %floor or ciel?

% TCSPC tail lenght in ns
% tcspc_tail_l = (tail_end - tail_start) * tcspc_bin_l;
tcspc_tail_l = (lt_end - lt_start) * lt_bin_l;

% Max Lifetime in ns
max_lt = 20;

%% SOFI Params

%duration of one frames (image) in seconds for SOFI JE 100 mu sec
img_d = 1e-3;

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


%number of events
n_events = numel(im_posx);

% measurment time in sec
m_time = head.ImgHdr_PixNum * head.ImgHdr_DwellTime;

%% ISM reasigment
fprintf('ISM reasigment\n');
%shift vectors from file negativ sign is already included
sv_file = matfile('shift_vector.m');
sv = sv_file.shift_vector;
shift_x     = sv(1, im_pix+1).';
shift_y     = sv(2, im_pix+1).';

sofi_sv = round(sv);

%apply ISM reassigment vektor
ISM_posx    = im_posx + shift_x;
ISM_posy    = im_posy + shift_y;

%counts per detector pixel
[pixel_int, ~]   = histcounts(im_pix,0:23);

%% Confocal Image
%generate confocal image
[sum_img, sum_lin, sum_size] = img_ps(im_posx, im_posy, s_pixl_x, s_pixl_y,1);

%remoce dark count
sum_img     = max(sum_img - sum(dc) * head.ImgHdr_PixelTime, 0);

%plot and save
reso_line_conf = [[105 105]; [30 65]];
img_plot(sum_img, spectrum, colf_name, 'confocal', 1, IM_R, reso_line_conf, 0, 1);

%% ISM Image
%crate ISM image
[ISM_img, ISM_lin, ISM_size]  = img_ps(ISM_posx, ISM_posy, s_pixl_x, s_pixl_y, ISM_binning);

%calculate ISM darkcount: Dc is linear and every ISM pixel gets count from
%23 pixels, either the same pixel in the sampel or a shifted on. But always
%23. Thus darkcount = sum(dc)
ISM_img  = max(ISM_img - sum(dc) * head.ImgHdr_PixelTime, 0);

reso_line_ISM = [[107 107]; [30 65]];
img_plot(ISM_img, spectrum, ISM_name, 'ISM', 1, IM_R, reso_line_ISM, 0, 1);

%% manuel lucy Richerson decon

if(deconv)
    PSF_file = matfile('PSF.m');
    psf = PSF_file.im;
    psf = psf/sum(psf, 'all');
    
    EID_file = matfile('EID.m');
    eid = EID_file.im;
    eid = eid/sum(eid, 'all');
    
    manuel_decon(ISM_img, eid, ...
                c_map = spectrum, s_name = ISM_docn_name, ...
                t_name = 'ISM + deconvolution', ...
                sb_lenght = 1,IM_R = IM_R,...
                reso = 0, save = 1);
end
%% Lifetime Image

if (lifetime)

    interval        = round(1.1* max(max(ISM_img)));
    n_pixel_ISM     = prod(ISM_size);
    ISM_lt          = zeros(n_pixel_ISM, lt_end-lt_start+1);

    h = waitbar(0,'binning');
    %find all photons in one ISM pixel, by seaching in an interval of max count
    %lengh + 1 
    lower = 1;
    upper = interval;

    % care y x könnten vertauscht sein
    % order photons in ISM lin index
    [ISM_lin, sort_index] = sort(ISM_lin);
    im_tcspc    = im_tcspc(sort_index);

    save_lower = zeros(n_pixel_ISM,1);
    save_upper = zeros(n_pixel_ISM,1);

    for i = 1:n_pixel_ISM
       [x,y] = ind2sub(ISM_size,i);
       ind = find(ISM_lin(lower:upper) == i);
       save_lower(i) = lower;
       save_upper(i) = upper;
       if ~isempty(ind)
           if(ISM_img(y,x) > lt_cut_off)
               [count, ~]   = histcounts(im_tcspc(lower + ind-1), 1:bin_factor:max_bin+1);
               count = count(lt_start:lt_end);
               ISM_lt(i,:)  = count./sum(count);
           end
%          set new upper edge to lower plus interval length
           upper    = min(lower + ind(end) + interval, n_events);
%          set lower edge to last found puls 1
           lower    = lower + ind(end);
       else
%            set new searche boundarys. because nothing was found keep
%            lower extend upper
           upper    = min(upper + interval, n_events);
       end
       if (mod(i,n_pixel_ISM/100) == 0)
            waitbar(i/n_pixel_ISM,h)
       end
    end
    close(h);
    fprintf('lifetime fit \n');

    %tail fit via pattern matching
    [ISM_lt_img,~]  = lt_patternMatching(ISM_lt, lt_bin_l, tcspc_tail_l, max_lt);

    %plot and save image
    ISM_lt_img  = reshape(ISM_lt_img, ISM_size);
%     ISM_lt_img = ISM_lt_img(150:190,130:180);
%     ISM_img = ISM_img(150:190,130:180);
    lt_img_plot(ISM_lt_img.', ISM_img, c_greenred, lt_cut_off, [.5 max_lt], 'Lifetime', LT_name, 2, IM_R, 1)
    
    figure
    histogram(ISM_lt_img(ISM_lt_img>0.01 & ISM_lt_img < 20),linspace(0.01,max_lt,100),'Normalization','probability')

    figure
    hold on
    plot(ISM_lin, 1:n_events, 'b-')
    plot( save_lower, 'r-')
    plot( save_upper, 'g-')
end
%% SOFI
if (sofi)

    interval        = round(1.1* max(max(sum_img)));
    n_pixel_sum     = prod(sum_size);
    SOFI_img        = zeros(n_pixel_sum, 1);
    SOFI_ism        = zeros(n_pixel_sum, 1);
    frame_times     = linspace(0, IM_dwell, n_frames + 1);

%     svLinInd        = sub2ind(sum_size,round(sv(2,:)),round(sv(1,:)));

    h = waitbar(0,'sofiing ?');
    %find all photons in one ISM pixel, by seaching in an interval of max count
    lower = 1;
    upper = interval;

    % care y x könnten vertauscht sein
    [sum_lin, sort_index] = sort(sum_lin);

    sofi_tcspc  = im_tcspc(sort_index);
    sofi_time   = im_time(sort_index);
    sofi_chan   = im_chan(sort_index);

    save_lower  = zeros(n_pixel_sum,1);
    save_upper  = zeros(n_pixel_sum,1);

    for i = 1:n_pixel_sum
       [y,x]            = ind2sub(sum_size,i);
       ind              = find(sum_lin(lower:upper) == i);
       save_lower(i)    = lower;
       save_upper(i)    = upper;
       if ~isempty(ind)

           if(1)%x == 23 && y == 50)
%                -1 ?
%               exact arrival time and detector channel of photons in
%               current image pixel(frame)
               f_time = sofi_time(lower + ind - 1) + sofi_tcspc(lower + ind-1) * time_R;
               f_chan = sofi_chan(lower + ind - 1);
               
               f_time = f_time - min(f_time);

%                t_max = max(f_time);
%                t_min = min(f_time);
% 
%                frame_time = t_min:img_d:t_max;
%                frames = size(frame_time,2)-1;

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
           end

           %set lower edge to last found puls 1
           n_lower      = lower + ind(end);
           %set new upper edge to lower plus interval length
           upper        = min(n_lower + interval, n_events);
       else
           %set new searche boundarys. because nothing was found keep lower
           %extend upper
           n_lower      = lower;
           upper        = min(upper + interval, n_events);
       end
       %set lower
       lower        = n_lower;

       if (mod(i,n_pixel_sum/100) == 0)
            waitbar(i/n_pixel_sum,h)
       end
    end
    close(h);
    fprintf('doing sofi \n');

    %plot and save image
    SOFI_img  = reshape(SOFI_img, sum_size);
    SOFI_ism  = reshape(SOFI_ism, sum_size);

    img_plot(SOFI_img, spectrum, 'sofi', 'sofi', 1, IM_R, reso_line_conf, 0, 0);
    img_plot(SOFI_ism, spectrum, 'sofi ism', 'sofi ism', 1, IM_R, reso_line_conf, 0, 0);

    figure
    hold on
    plot(sum_lin, 1:n_events, 'b-')
    plot( save_lower, 'r-')
    plot( save_upper, 'g-')
end

%% Additional figures for controle

% Normalised monoexponetial decay with background
pfun_monoexp = @(tau,x) exp(-x(:)./tau); 
pfun_monoexpBG = @(tau,b,x)b./numel(x(:))+(1-b).*pfun_monoexp(tau,x)./sum(pfun_monoexp(tau,x),1);

% figure
% [c_count, ~] = histcounts(im_tcspc, 1:max_bin+1);
% plot(c_count(1:tail_end))
% set(gca, 'YScale', 'log')

% hex_plot(pixel_int.',hot)
% plot_pixeldecay(im_tcspc,im_chan);


% dt      = tcspc_bin_l;
% p       = max_bin * tcspc_bin_l;
% [irf,~] = histcounts(irf_tcspc,1:max_bin+1);
% [y, ~]  = histcounts(im_tcspc,1:max_bin+1);
% irf     = irf(1:tail_end);
% % y       = max(y - sum(dc) * head.MeasDesc_AcquisitionTime * 1.6e-3/max_bin,1e-1);
% y_all       = y(1:tail_end);
% max_y   = max(y_all);
% max_irf = max(irf);
% irf     = irf.*max_y./max_irf;
% taus    = [5, 12];
% lim     = [0 0; 25 25];

% [count, ~]   = histcounts(im_tcspc, 1:bin_factor:max_bin+1);
% y_tail =  count(lt_start:lt_end);
% y_tail = y_tail./sum(y_tail);
% [over_all_lt,over_all_back]  = lt_patternMatching(y_tail, lt_bin_l, tcspc_tail_l, max_lt);
% disp(over_all_lt)
% 
% xx = 0:(lt_end-lt_start);
% tcspc_t = 0:lt_bin_l:tcspc_tail_l;
% yy = pfun_monoexpBG(over_all_lt,over_all_back,tcspc_t);
% figure
% hold on
% % [c_count, ~] = histcounts(y, 1:max_bin+1);
% plot(y_tail)
% plot(yy.')
% legend('data','pattern Matching')
% set(gca, 'YScale', 'log')

% [c, offset, A, tau, ~, ~, ~, ~, ~, ~] = Fluofit(irf, y_all, p, dt, taus, lim, 1,1);

