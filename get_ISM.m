%% Clear shit  up
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F)
clear
close all
clc

%%

lifetime    = 0;
s_lifetime  = 0;
d_lifetime  = 0;
t_lifetime  = 0;
q_lifetime  = 0;
deconv      = 0;
sofi        = 0;
add_plt     = 0;
%% Loadind data
%rgb values for color map black,blue,cyan,green,yellow,orange?,red,magenta
sp1     = 1:255/7:256;
sp2     = [[0 0 0]; [0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]; [1 0 1]];
br1     = 1:255/5:256;
br2     = [[0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]];
grl     = 1:255/3:256;
gr2     = [[0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]];
gl      = 1:255/1:256;
% g2      = [[0 0 0]; [0.454325 0.806075 0.332481]];
g2      = [[0 0 0];[0 1 0]];
rl      = 1:255/1:256;
r2      = [[0 0 0]; [1.000000 0.6 0.5]];
% r2      = [[0 0 0]; [1 0 0]];
bl      = 1:255/1:256;
b2      = [[0 0 0]; [0.000000 0.8 1.000000]];
% b2      = [[0 0 0]; [0 0 1]];
y2      = [[0 0 0]; [1 1 0]];

lambda  = 1:256;

spectrum    = interp1(sp1, sp2, lambda);
c_greenred  = interp1(grl, gr2, lambda);
c_bluered   = interp1(br1, br2, lambda);
c_green     = interp1(gl, g2, lambda);
c_red       = interp1(rl, r2, lambda);
c_blue      = interp1(bl, b2, lambda);
c_yellow    = interp1(bl, y2, lambda);
c_map       = cmap_isoluminant75;

fname   = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230207\neurons_g4_af647_gm130_one_pmt_006.ptu';
dcname  = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230203\dc_m15_pmt_001.ptu';
irfname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220816\irf_ex470nm_005.ptu';

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
[im_chan,im_tcspc,im_posy,im_posx,im_time,head] = read_ISM(fname);

ind             = im_posy<=head.ImgHdr_PixY;
im_tcspc        = im_tcspc(ind);
im_chan         = im_chan(ind);
im_posx         = im_posx(ind);
im_posy         = im_posy(ind);
im_time         = im_time(ind);
im_posx         = im_posx * head.ImgHdr_PixX;
im_time         = im_time./head.TTResult_SyncRate; % photon arrival in seconds

%remove time delay due to unsyncronised Multi Harps
im_tcspc = remove_MHH_offset(im_tcspc,im_chan, head.max_bin, 500);

%Dark count in count per second
[dc, bin_dc] = get_DC(dcname,1);

%get IRF to set tail for tailfit
% [fwhm, indMax, ~, ~,irf_tcspc] = get_IRF(irfname, 1);
%% Parameters and magic numbers(please fix) AND FIX NAMEN FÜR TCSPC

%number of pixels in recorded image
s_pixl_x    = head.ImgHdr_PixX;
s_pixl_y    = head.ImgHdr_PixY;

%scan resolution step per pixel in µm
IM_R = head.ImgHdr_PixResol;

%ISM Binning UGO[1.02] PQ[1.04] UGO220302[1.01] 20nm[1.0075] g4mix
%[1.01/1.0075]
ISM_binning = 1.01;

%temporal resolution of TCSPC in seconds
time_R = mean(head.MeasDesc_Resolution);

%number of pixels of the detector
if(sum(head.HWInpChan_Enabled) > 23)
    n_pixl = 32;
else
    n_pixl = 23;
end

%number of events
n_events = numel(im_chan);

max_bin = head.max_bin;

% TCSPC binning factor
bin_factor = 20;

% minimum numbers of photons to calculate lifetime
lt_cut_off = 100;

% tcspc tail_start in tcspc bin number !!FIX NEEDED!!
% tail_start = indMax + fwhm; diode [600/700,275], TiSa [250], SEPIA [1700]
tail_start = 700;

% tcspc tail_end in tcspc bin number
tail_end = max_bin - 100;

[tail_t, tail_bin_l, tail_start_time, tcspc_t] = get_lifetime_bins(tail_start, tail_end, time_R, bin_factor);

% Max Lifetime in ns
max_lt = 8;

% Normalised monoexponetial decay with background
pfun_monoexp = @(tau,x,dt) dt(:).*exp(-x(:)./tau); 
pfun_monoexpBG = @(tau,b,x,dt)b./numel(x(:))+(1-b).*pfun_monoexp(tau,x,dt)./sum(pfun_monoexp(tau,x,dt),1);

pixel_img   = histcounts(im_chan,0:n_pixl);
fiber_img   = zeros(37,1);
fiber_img(fiber_code()+1)   = pixel_img;
hex_plot(fiber_img,hot)

%% Multiple tau scale

m_tau = 0;
if(m_tau)
    % bin size
    tau_s = @(delta, i) delta*2.^(floor(i/8));
    
    % number of bin requierd for given tau scale. approx sum 2^(floor(i/8))
    % with sum 2^(i/8-1/2) and use geometric sum
    
    % constant form geometric series
    geometric_const = 2^(3/8)*(nthroot(2,8) - 1);
    N = floor(8 * log2(tail_l/bin_factor * geometric_const));
    
    % tcspc bins according to multi tau scale
    tcspc_bin   = tau_s(bin_factor, 1:N);
    % and bin edges
    tcspc_t     = tail_start + cumsum([0, tcspc_bin(1:end-1)]);

    % and in ns for tail fit
    tail_bin_l    = tcspc_bin_l*tcspc_bin(1:end-1);
    tail_t      = tcspc_bin_l*tcspc_t(1:end-1);

%else
    tcspc_bin   = bin_factor;
    tcspc_t     = tail_start:bin_factor:tail_end;

    % and in ns for tail fit
    tail_bin_l    = tcspc_bin_l*tcspc_bin;
    tail_t      = tcspc_bin_l*tcspc_t(1:end-1);
end

%% Dark Count

% bin_dc = sum(bin_dc, 1);
% 
% bin_dc = movsum(bin_dc, [0 bin_factor-1]);
% 
% bin_dc = bin_dc(tcspc_t);

%% ISM reasigment
fprintf('ISM reasigment\n');
%shift vectors from file negativ sign is already included
sv_file = matfile('shift_vector.m');
sv = sv_file.shift_vector;
shift_x     = sv(1, im_chan+1).';
shift_y     = sv(2, im_chan+1).';

%apply ISM reassigment vektor
ISM_posx    = im_posx + shift_x./(IM_R/0.05);
ISM_posy    = im_posy + shift_y./(IM_R/0.05);

clear shift_y shift_x;

%% Confocal Image
%generate confocal image
[sum_img, sum_lin, sum_size] = img_ps(im_posx, im_posy, s_pixl_x, s_pixl_y,1);

% sum_img = sum_img(154:204,260:312);

clear im_posy im_posx;
%remoce dark count
%sum_img = max(sum_img - sum(dc) * head.ImgHdr_PixelTime, 0);

%plot and save
reso_line_conf = [[42 92]; [98 98]];
img_plot(max(sum_img - sum(dc) * head.ImgHdr_PixelTime, 0), hot, colf_name, 'confocal', 1, IM_R, 1, reso_line_conf, 0, 1);

%% ISM Image
%crate ISM image
[ISM_img, ISM_lin, ISM_size]  = img_ps(ISM_posx, ISM_posy, s_pixl_x, s_pixl_y, ISM_binning);

% ISM_img = ISM_img(156:206,262:312);

% clear ISM_posy ISM_posx;
%calculate ISM darkcount: Dc is linear and every ISM pixel gets count from
%23 pixels, either the same pixel in the sampel or a shifted on. But always
%23. Thus darkcount = sum(dc)

%plot and save
reso_line_ISM = [[42 92]; [98 98]];
img_plot(max(ISM_img - sum(dc) * head.ImgHdr_PixelTime, 0), hot, ISM_name, 'ISM', 1, IM_R, ISM_binning, reso_line_ISM, 0, 1);

%% manuel lucy Richerson decon
if(deconv)
    PSF_file = matfile('PSF.m');
    psf      = PSF_file.im;
    psf      = psf/sum(psf, 'all');
    
    EID_file = matfile('EID.m');
    eid      = EID_file.im;
    eid      = eid/sum(eid, 'all');
    
    manuel_decon(max(ISM_img - sum(dc) * head.ImgHdr_PixelTime, 0), eid, ...
                c_map = hot, s_name = ISM_docn_name, ...
                t_name = 'ISM + deconvolution', ...
                sb_lenght = 1, IM_R = IM_R, pix_bin = ISM_binning, ...
                reso = 1, save = 1);
end

%% Lifetime Image
if (lifetime)

    interval        = ceil(1.01* max(max(ISM_img)));
    n_pixel_ISM     = prod(ISM_size);
    ISM_lt          = zeros(n_pixel_ISM, numel(tail_t));
    ISM_lt_amp      = zeros(n_pixel_ISM, 4);
    ISM_lt_rgb      = zeros([ISM_size 3]);

    % generate decay patterns from FL selecion[1.37 2.36 3.0]
    pattern_tau = [1.4 2.4 3.4 inf];
    pattern = pfun_monoexpBG(pattern_tau, 0, tail_t-tail_start_time, tail_bin_l);

    ISM_lt_short_name   = compose('%s lt short %0.1f ns', ISM_name, pattern_tau(1)); %1.37
    ISM_lt_middle_name  = compose('%s lt middle %0.1f ns', ISM_name, pattern_tau(2)); %2.37
    ISM_lt_long_name    = compose('%s lt long %0.1f ns', ISM_name, pattern_tau(3)); % 3.0

    h = waitbar(0,'binning');
    %find all photons in one ISM pixel, by seaching in an interval of max count
    %lengh + 1 
    lower = 1;
    upper = interval;
    
    % Sort for faster seache
    [ISM_lin, sort_index] = sort(ISM_lin);
    im_tcspc    = im_tcspc(sort_index);

    for i = 1:n_pixel_ISM
       [y,x]    = ind2sub(ISM_size,i);
       ind      = find(ISM_lin(lower:upper) == i);
       if ~isempty(ind)
           if(ISM_img(y,x) >= lt_cut_off)
               [count, ~]   = histcounts(im_tcspc(lower + ind-1), tcspc_t);
               ISM_lt_amp(i,:) = lsqnonneg(pattern, count.');
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
       if (mod(i,n_pixel_ISM/100) < 1)
            waitbar(i/n_pixel_ISM,h)
       end
    end
    close(h);
    fprintf('lifetime fit \n');

    ISM_lt_amp_img = reshape(ISM_lt_amp, [ISM_size(1), ISM_size(2), 4]);
    img_plot(ISM_lt_amp_img(:,:,1), c_blue, ISM_lt_short_name{1}, 'Cy2: SYT 1', 4, IM_R, ISM_binning, reso_line_ISM, 0, 1);
    img_plot(ISM_lt_amp_img(:,:,2), c_red, ISM_lt_middle_name{1}, 'Alexa: TOMM20', 4, IM_R, ISM_binning, reso_line_ISM, 0, 1);
    img_plot(ISM_lt_amp_img(:,:,3), c_green, ISM_lt_long_name{1}, 'OG: GFAP', 4, IM_R, ISM_binning, reso_line_ISM, 0, 1);
    
    % adjust contrast
    P = prctile(ISM_lt_amp_img(ISM_lt_amp_img(:,:,2)>0),[10, 95], 'all');
    ISM_lt_rgb(:,:,1) = adapthisteq(ISM_lt_amp_img(:,:,2)./max(ISM_lt_amp_img(:,:,2),[],'all'),'ClipLimit',0.005);%red

    P = prctile(ISM_lt_amp_img(ISM_lt_amp_img(:,:,3)>0),[10, 99.7], 'all');
    ISM_lt_rgb(:,:,2) = adapthisteq(ISM_lt_amp_img(:,:,3)./max(ISM_lt_amp_img(:,:,3),[],'all'),'ClipLimit',0.005);%green

    P = prctile(ISM_lt_amp_img(ISM_lt_amp_img(:,:,1)>0),[10, 99.7], 'all');
    ISM_lt_rgb(:,:,3) = adapthisteq(ISM_lt_amp_img(:,:,1)./max(ISM_lt_amp_img(:,:,1),[],'all'),'ClipLimit',0.005);%blue

    img_plot(ISM_lt_rgb(:,:,3), c_blue, ISM_lt_short_name{1}, 'Cy2: SYT 1', 4, IM_R, ISM_binning, reso_line_ISM, 0, 0);
    img_plot(ISM_lt_rgb(:,:,1), c_red, ISM_lt_middle_name{1}, 'Alexa: TOMM20', 4, IM_R, ISM_binning, reso_line_ISM, 0, 0);
    img_plot(ISM_lt_rgb(:,:,2), c_green, ISM_lt_long_name{1}, 'OG: GFAP', 4, IM_R, ISM_binning, reso_line_ISM, 0, 0);

    comb_name = compose('%s multicolor %0.1f %0.1f %.01f', LT_name, pattern_tau(1), pattern_tau(2), pattern_tau(3));

    rgb_img_plot(ISM_lt_rgb, comb_name{1}, 'Red: PSD95, Green: GFAP, Blue: SYT 1', 4, IM_R, 1)

%     test_img = ISM_lt_amp_img(:,:,1:3);
%     test_img = ISM_lt_amp_img./repmat(max(ISM_lt_amp_img,[],[1 2]), [size(ISM_lt_amp_img,[1 2]) 1]);
%     test_img = test_img(:,:,1:3);
%     test_color =[g2(2,:); r2(2,:); b2(2,:)];
%     test_color =[[1 0 1]; [0 1 1]; [0.96 1 0]];
% 
%     test_rgb = tensorprod(test_img,test_color,3,1);
%     test_rgb = test_rgb./repmat(max(test_rgb,[],[1 2]), [size(test_rgb,[1 2]) 1]);
%     rgb_img_plot(test_rgb, LT_name, 'Red: GFAP, Green: SYT, Blue: PSD95', 4, IM_R, 0)
end

%% Single Expoential lifetime
if(s_lifetime)
    %get pixel lifetimes via MLE pattern matching
    lt_img = get_single_lifetime(ISM_img,ISM_lin,...
        im_tcspc = im_tcspc, tail_t = tail_t, tail_bin_l = tail_bin_l,...
        tail_start_time = tail_start_time, tcspc_t = tcspc_t,...
        max_lt = max_lt, lt_cut_off = lt_cut_off, fname = LT_name);

    %plot pxel wise lt values scaled with image intensity
    lt_img_plot(lt_img, ISM_img, c_map, lt_cut_off, [0.2 5], 'Lifetime', LT_name, 4, IM_R, 1)   
end

%% Tripple lifetime unmixing
if(t_lifetime)
    % generate decay patterns from FL selecion[1.37 2.36 3.0]
    pattern_tau = [0.34 1.16 3.38 inf];
    pattern = pfun_monoexpBG(pattern_tau,0,tail_t-tail_start_time,tail_bin_l);

    [lt_amp_img] = get_unmixed_lifetime(img,lin,im_tcspc = im_tcspc, ...
        tail_t = tail_t, tail_bin = tail_bin_l, tail_start_time = tail_start_time, ...
        tcspc_t = tcspc_t, pattern = pattern, name = ISM_name);

    %plot single images

    %function that optimises single images and pots merged image
end

%% Double lifetime unmixing
if(d_lifetime)
    % generate decay patterns from FL selecion[1.37 2.36 3.0]
    pattern_tau = [1.2 2.7 inf];
    pattern = pfun_monoexpBG(pattern_tau,0,tail_t-tail_start_time,tail_bin_l);

    % Struckture and Fluorophore name
    names = {'SYT1', 'GFAP'; 'Cy2', 'Oregon Green 488'};

    lt_short_name   = compose('%s lt short %0.1f ns', ISM_name, pattern_tau(1)); %1.37
    lt_long_name    = compose('%s lt long %0.1f ns', ISM_name, pattern_tau(2)); % 3.0

    shot_title = append(names{1,1}, ': ', names{2,1});
    long_title = append(names{1,2}, ': ', names{2,2});
    
    %unmix the two patterns
    [lt_amp_img] = get_bi_lifetime(ISM_img,ISM_lin,im_tcspc = im_tcspc, ...
        tail_t = tail_t, tail_bin = tail_bin_l, tail_start_time = tail_start_time, ...
        tcspc_t = tcspc_t, pattern = pattern, pattern_tau = pattern_tau, ...
        name = names, fname = ISM_name);
    
    %plot individual strucktures
    img_plot(lt_amp_img(:,:,1), c_blue, lt_short_name{1}, shot_title, 4, IM_R, ISM_binning, reso_line_ISM, 0, 1);
    img_plot(lt_amp_img(:,:,2), c_green, lt_long_name{1}, long_title, 4, IM_R, ISM_binning, reso_line_ISM, 0, 1);

    %merge color plots
    merge_two_lt(lt_amp_img, names, LT_name, pattern_tau, IM_R);

end

%% Double lifetime double color
if(q_lifetime)
    %start bin of second color
    red_start = 5000;

    blue = im_tcspc < red_start;
    red  = im_tcspc >= red_start;

    blue_tcspc      = im_tcspc(blue);
    blue_chan       = im_tcspc(blue);
    blue_posx       = im_tcspc(blue);
    blue_posy       = im_tcspc(blue);
    blue_ISM_posx   = ISM_posx(blue);
    blue_ISM_posy   = ISM_posy(blue);

    red_tcspc      = im_tcspc(red);
    red_chan       = im_tcspc(red);
    red_posx       = im_tcspc(red);
    red_posy       = im_tcspc(red);
    red_ISM_posx   = ISM_posx(red);
    red_ISM_posy   = ISM_posy(red);

    %crate blue ISM image
    [blue_ISM_img, blue_ISM_lin, ~]  = img_ps(blue_ISM_posx, blue_ISM_posy, s_pixl_x, s_pixl_y, ISM_binning);
    
    %plot and save
    reso_line_ISM = [[107 107]; [30 65]];
    img_plot(blue_ISM_img, hot, ISM_name, 'blue ISM', 4, IM_R, ISM_binning, reso_line_ISM, 0, 0);

    %crate blue ISM image
    [red_ISM_img, red_ISM_lin, ~]  = img_ps(red_ISM_posx, red_ISM_posy, s_pixl_x, s_pixl_y, ISM_binning);
    
    %plot and save
    reso_line_ISM = [[107 107]; [30 65]];
    img_plot(red_ISM_img, hot, ISM_name, 'red ISM', 4, IM_R, ISM_binning, reso_line_ISM, 0, 0);

    %lifetiem analysis
    [tail_t, tail_bin_l, tail_start_time, tcspc_t] = get_lifetime_bins(335, 4982, time_R, bin_factor);
    blue_pattern_tau = [1.5 3.4 inf];
    blue_pattern = pfun_monoexpBG(blue_pattern_tau,0,tail_t-tail_start_time,tail_bin_l);

    % Struckture and Fluorophore name
    names = {'TOMM20', 'GFAP'; 'Cy2', 'OG 488'};

    lt_short_name   = compose('%s blue lt short %0.1f ns', ISM_name, blue_pattern_tau(1)); %1.37
    lt_long_name    = compose('%s blue lt long %0.1f ns', ISM_name, blue_pattern_tau(2)); % 3.0

    shot_title = append(names{1,1}, ': ', names{2,1});
    long_title = append(names{1,2}, ': ', names{2,2});
    
    %unmix the two patterns
    [lt_amp_img] = get_bi_lifetime(blue_ISM_img, blue_ISM_lin, im_tcspc = blue_tcspc, ...
        tail_t = tail_t, tail_bin = tail_bin_l, tail_start_time = tail_start_time, ...
        tcspc_t = tcspc_t, pattern = blue_pattern, pattern_tau = blue_pattern_tau, ...
        name = names, fname = ISM_name);
    
    %plot individual strucktures
    img_plot(lt_amp_img(:,:,1), c_blue, lt_short_name{1}, shot_title, 4, IM_R, ISM_binning, reso_line_ISM, 0, 1);
    img_plot(lt_amp_img(:,:,2), c_green, lt_long_name{1}, long_title, 4, IM_R, ISM_binning, reso_line_ISM, 0, 1);

% %     get pixel lifetimes via MLE pattern matching
%     lt_img = get_single_lifetime(blue_ISM_img, blue_ISM_lin, ...
%         im_tcspc = blue_tcspc, tail_t = tail_t, tail_bin_l = tail_bin_l,...
%         tail_start_time = tail_start_time, tcspc_t = tcspc_t,...
%         max_lt = max_lt, lt_cut_off = lt_cut_off);
% 
%     %plot pxel wise lt values scaled with image intensity
%     lt_img_plot(lt_img, blue_ISM_img, c_map, lt_cut_off, [0.2 4], 'Lifetime', LT_name, 4, IM_R, 1) 


    %lifetiem analysis red
    [tail_t, tail_bin_l, tail_start_time, tcspc_t] = get_lifetime_bins(5373, 9990, time_R, bin_factor);
    red_pattern_tau = [1.7 3.8 inf];
    red_pattern = pfun_monoexpBG(red_pattern_tau,0,tail_t-tail_start_time,tail_bin_l);

    % Struckture and Fluorophore name
    names = {'GM130', 'PSD95'; 'AF 647', 'AS 635P'};

    lt_short_name   = compose('%s red lt short %0.1f ns', ISM_name, red_pattern_tau(1)); %1.37
    lt_long_name    = compose('%s red lt long %0.1f ns', ISM_name, red_pattern_tau(2)); % 3.0

    shot_title = append(names{1,1}, ': ', names{2,1});
    long_title = append(names{1,2}, ': ', names{2,2});
    
    %unmix the two patterns
    [lt_amp_img] = get_bi_lifetime(red_ISM_img, red_ISM_lin, im_tcspc = red_tcspc, ...
        tail_t = tail_t, tail_bin = tail_bin_l, tail_start_time = tail_start_time, ...
        tcspc_t = tcspc_t, pattern = red_pattern, pattern_tau = red_pattern_tau, ...
        name = names, fname = ISM_name);
    
    %plot individual strucktures
    img_plot(lt_amp_img(:,:,1), c_yellow, lt_short_name{1}, shot_title, 4, IM_R, ISM_binning, reso_line_ISM, 0, 1);
    img_plot(lt_amp_img(:,:,2), c_red, lt_long_name{1}, long_title, 4, IM_R, ISM_binning, reso_line_ISM, 0, 1);

%     get pixel lifetimes via MLE pattern matching
%     lt_img = get_single_lifetime(red_ISM_img, red_ISM_lin, ...
%             im_tcspc = red_tcspc, tail_t = tail_t, tail_bin_l = tail_bin_l,...
%             tail_start_time = tail_start_time, tcspc_t = tcspc_t,...
%             max_lt = max_lt, lt_cut_off = lt_cut_off, fname = LT_name);

    %plot pxel wise lt values scaled with image intensity
%     lt_img_plot(lt_img, red_ISM_img, c_map, lt_cut_off, [0.2 4], 'Lifetime', LT_name, 4, IM_R, 1) 
end

%% SOFI
if (sofi)

    ISM_SOFI_img = get_ISMSOFI(ISM_img, ISM_lin, im_time, head);

    img_plot(ISM_SOFI_img, spectrum, 'sofi ism', 'sofi ism', 1, IM_R, ISM_binning, reso_line_conf, 0, 1);

    SOFI_img = get_ISMSOFI(sum_img, sum_lin, im_time, head);

    img_plot(SOFI_img, spectrum, 'sofi sum', 'sofi sum', 1, IM_R, 1, reso_line_conf, 0, 1);
end

%% Additional figures for controle
if(add_plt)    
    f = figure;
    ax = axes(f);
    histogram(ax, im_tcspc, 1:max_bin, 'EdgeAlpha',0)
    set(ax,'YScale', 'log')
    xlabel('tcspc bin')
    ylabel('count')
    
    [count, edges]   = histcounts(im_tcspc, tcspc_t);
    
    [tripple_amp,tau] = tripple_exp(tail_t.'-tail_start_time,count.',0);
    disp('over all tau and amp')
    disp(tau)
    disp(tripple_amp./sum(tripple_amp))
    
    xx = tail_t-tail_start_time;
    count = count./sum(count);
    [over_all_lt,over_all_back]  = lt_patternMatching(count, xx, tail_bin_l, max_lt);
    disp(over_all_lt)
    
    h = figure;
    ax = axes(h);
    
    ax1 = subplot(3,4,[1,2,3,5,6,7,9,10,11]);
    p1 = plot(xx, count, 'o', 'MarkerSize', 10, 'Color', '#EDB120');
    hold on
    
    yy = pfun_monoexpBG(over_all_lt, over_all_back, tail_t, tail_bin_l);
    p2 = plot(tail_t-tail_start_time, yy.','g-','LineWidth',2);
    
    amp1 = tripple_amp(1);
    amp2 = tripple_amp(2);
    amp3 = tripple_amp(3);
    offset = tripple_amp(4);
    y4 =  amp1*exp(-xx/tau(1)) + amp2*exp(-xx/tau(2)) + amp3*exp(-xx/tau(3)) + offset;
    y4 = y4./sum(y4);
    p3 = plot(xx,y4,'-','LineWidth',2, 'Color', '#A2142F');
    
    amp1 = sum(ISM_lt_amp_img(:,:,1), 'all');
    amp2 = sum(ISM_lt_amp_img(:,:,2), 'all');
    amp3 = sum(ISM_lt_amp_img(:,:,3), 'all');
    offset = sum(ISM_lt_amp_img(:,:,4), 'all');
    y3 =  amp1*pattern(:,1) + amp2*pattern(:,2) + amp3*pattern(:,3) + offset*pattern(:,4);
    y3 = y3./sum(y3);
    plot(xx,y3,'b-','LineWidth',2)
    
    legend('tail decay','pattern Matching', 'true tripple', 'from image tripple')
    % legend('tail decay','tripple exponential fit')
    set(ax1, 'YScale', 'log')
    % ylim([0.5*min(yy) inf])
    xlim('padded')
    ylim('padded')
    xlabel('time [ns]')
    ylabel('normalized count')
    title('TCSPC decay of an individual pixel');
    set(findall(h,'-property','FontSize'),'FontSize',15)
    
    ax2 = subplot(3,4,4);
    plot(tail_t-tail_start_time, pattern(:,1)./max(pattern(:,1)),'b-')
    set(ax2, 'YScale', 'log')
    ylim([1e-4 1])
    xlim([0 20])
    title('short');
    
    ax3 = subplot(3,4,8);
    plot(tail_t-tail_start_time, pattern(:,2)./max(pattern(:,2)),'r-')
    set(ax3, 'YScale', 'log')
    ylim([1e-4 1])
    xlim([0 20])
    title('middle');
    
    ax4 = subplot(3,4,12);
    plot(tail_t-tail_start_time, pattern(:,3)./max(pattern(:,3)),'g-')
    set(ax4, 'YScale', 'log')
    ylim([1e-4 1])
    xlim([0 20])
    title('long');
    
    disp('amps from pixel fit')
    fprintf('%f %f %f %f\n',amp1/(amp1+amp2+amp3+offset), amp2/(amp1+amp2+amp3+offset),...
        amp3/(amp1+amp2+amp3+offset), offset/(amp1+amp2+amp3+offset))
    
    plot_pixeldecay(im_tcspc,im_chan,32)

    % dt      = tcspc_bin_l;
    % p       = max_bin * tcspc_bin_l;
    % [irf,~] = histcounts(irf_tcspc,1:max_bin+1);
    % [y, ~]  = histcounts(im_tcspc,1:max_bin+1);
    % irf     = irf(1:tail_end);
    % % y     = max(y - sum(dc) * head.MeasDesc_AcquisitionTime * 1.6e-3/max_bin,1e-1);
    % y_all   = y(1:tail_end);
    % max_y   = max(y_all);
    % max_irf = max(irf);
    % % irf     = irf.*max_y./max_irf;
    % taus    = [0.5 1 2 3.45];
    % lim     = [0 0 0 0;8 8 8 8];
    % 
    % [c, offset, A, tau, ~, ~, ~, ~, ~, ~] = Fluofit(irf, y_all, p, dt, taus, lim, 0, 1);

end