function get_ISM(fname, dcname, irfname, options)
%% Clear shit  up
% F = findall(0,'type','figure','tag','TMWWaitbar');
% delete(F)
% clear
% close all
% clc

%%

arguments
    fname char             
    dcname char           
    irfname char          
    options (1,1) struct
%     options.lifetime     double = 0
%     options.s_lifetime double   = 0
%     options.d_lifetime double   = 0
%     options.d_color    double   = 0
%     options.t_lifetime double   = 0
%     options.q_lifetime double   = 0
%     options.deconv     double   = 0
%     options.frw        double   = 0
%     options.sofi       double   = 0
%     options.add_plt    double   = 0
%     options.save_image double   = 0
%     options.plot_reso  double   = 0
end
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

% fname   = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230413_CMs_fixed\CM_fixed_20kPa_ACTN2_citrine_Paxilin_SPAD_019.ptu';
% dcname  = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230314_CMs_life\dc_spad_001.ptu';
% irfname = 'D:\PHD\Data\2022\220816\irf_ex470nm_005.ptu';

%Image title and name from file name
tmp_name        = strsplit(fname, '\');
folder_name     = tmp_name{end-1};
date            = strsplit(folder_name, '_');
date            = date{1};
img_name        = tmp_name{end};
img_name        = strsplit(img_name, '.');
img_name        = img_name{end-1};
img_name        = append(date,' ',img_name);
ISM_name        = append(img_name, ' ISM');
ISM_fw_name     = append(img_name, ' ISM Fourier-reweighted');
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
% im_tcspc = remove_MHH_offset(im_tcspc,im_chan, head.max_bin, 500);

%Dark count in count per second
[dc, ~] = get_DC(dcname,0);

%get IRF to set tail for tailfit
% [fwhm, indMax, ~, ~,irf_tcspc] = get_IRF(irfname, 0);
%% Parameters and magic numbers(please fix) AND FIX NAMEN FÜR TCSPC

%number of pixels in recorded image
s_pixl_x    = head.ImgHdr_PixX;
s_pixl_y    = head.ImgHdr_PixY;

%scan resolution step per pixel in µm
IM_R = head.ImgHdr_PixResol;

%ISM Binning UGO[1.02] PQ[1.04] UGO220302[1.01] 20nm[1.0075] g4mix
%[1.01/1.0075]pmt
ISM_binning = 1.01;

%temporal resolution of TCSPC in seconds
time_R = mean(head.MeasDesc_Resolution);

%number of pixels of the detector
if(sum(head.HWInpChan_Enabled) > 23)
    n_pixl = 32;
else
    n_pixl = 23;
    %need fix for switch box
    im_chan = im_chan - 9;
end

max_bin = head.max_bin;

% TCSPC binning factor
bin_factor = 20;

% minimum numbers of photons to calculate lifetime
lt_cut_off = 10;

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

%% ISM reasigment
fprintf('ISM reasigment\n');
%shift vectors from file negativ sign is already included
if n_pixl == 23
    sv_file_name = 'SPAD_shift_vectors.m';
else
    sv_file_name = 'MPMT_shift_vectors.m';
end

sv_file = matfile(sv_file_name);
sv = sv_file.shift_vector;
shift_x     = sv(1, im_chan+1).';
shift_y     = sv(2, im_chan+1).';

%apply ISM reassigment vektor
ISM_posx    = im_posx + shift_x;%./(IM_R/0.05);
ISM_posy    = im_posy + shift_y;%./(IM_R/0.05);

clear shift_y shift_x;

%% Confocal Image
%generate confocal image
[sum_img, sum_lin, ~] = img_ps(im_posx, im_posy, s_pixl_x, s_pixl_y,1);

% ind_x = im_posx >= 726.5 & im_posx < 727.5;
% ind_y = im_posy == 154;
% 
% ind = ind_y&ind_x;
% 
% det_pixel = im_chan(ind);
% 
% [pix_count,~] = histcounts(det_pixel, n_pixl);
% hex_plot(pix_count, hot);
% 
% det_x = -1.*[-2,-1,0,1,2,-3/2,-1/2,1/2,3/2,-2,-1,0,1,2,-3/2,-1/2,1/2,3/2,-2,-1,0,1,2];
% A = ones(1,5).*sqrt(3);
% B = zeros(1,5);
% det_y = 1.*[A,A(1:end-1)./2,B,-A(1:end-1)./2,-A];
% 
% d_gauss = @(a, b, s1, s2, x, y) ...
%     a * exp(-(x.^2+y.^2)/(2*s1^2) + b * exp(-(x.^2+y.^2)/(2*s2^2)));
% d_gauss_fitt = fit([det_x.',det_y.'], pix_count.', d_gauss, ...
%     'StartPoint', [10, 10, 0.7, 1.7], ...
%     'Lower', [1, 1, 0.1, 0.1]);
% [xx, yy] = meshgrid(-2:0.1:2,-2:0.1:2);
% figure
% hold on
% test = d_gauss(d_gauss_fitt.a,d_gauss_fitt.b,d_gauss_fitt.s1,d_gauss_fitt.s2,xx,yy);
% scatter3(det_x,det_y, pix_count)
% surf(xx,yy,test)
% axis square
% clear im_posy im_posx;

%plot and save
reso_line_conf = [[110 110]; [29 79]];
sum_img_dc = max(sum_img - sum(dc) * head.ImgHdr_PixelTime, 0);
img_plot( sum_img_dc, hot, colf_name, ...
    'confocal', 1, IM_R, 1, reso_line_conf, options.plot_reso, options.save_image);

%% ISM Image
%crate ISM image
[ISM_img, ISM_lin, ~]  = img_ps(ISM_posx, ISM_posy, s_pixl_x, s_pixl_y, ISM_binning);

%calculate ISM darkcount: Dc is linear and every ISM pixel gets count from
%23 pixels, either the same pixel in the sampel or a shifted on. But always
%23. Thus darkcount = sum(dc)

%plot and save
reso_line_ISM = [[111 111]; [29 79]];
img_plot(max(ISM_img - sum(dc) * head.ImgHdr_PixelTime, 0), hot, ISM_name, ...
    'ISM', 1, IM_R, ISM_binning, reso_line_ISM, options.plot_reso, ...
    options.save_image);

%% Fourier-reweighted ISM

if options.frw
    %calculate Fourier-reweighted ISM image
    W_ISM_img1 = f_reweighting( max(ISM_img - sum(dc) * head.ImgHdr_PixelTime, 0), IM_R);

    %plot and save fr ISM
    reso_line_frw = [[111 111]; [29 79]];
    img_plot(W_ISM_img1, hot, ISM_fw_name, 'Fourier-reweighted ISM', ...
        1, IM_R, ISM_binning, reso_line_frw, options.plot_reso, ...
        options.save_image);
end

%% Single Expoential lifetime
if options.s_lifetime
    %get pixel lifetimes via MLE pattern matching
    lt_img = get_single_lifetime(ISM_img,ISM_lin,...
        im_tcspc = im_tcspc, tail_t = tail_t, tail_bin_l = tail_bin_l,...
        tail_start_time = tail_start_time, tcspc_t = tcspc_t,...
        max_lt = max_lt, lt_cut_off = lt_cut_off, fname = LT_name);

    %plot pxel wise lt values scaled with image intensity
    lt_img_plot(lt_img, ISM_img, c_map, lt_cut_off, [0.2 5], 'Lifetime', LT_name, 4, IM_R, 1)   
end

%% Tripple lifetime unmixing
if options.t_lifetime
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
if options.d_lifetime
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
if(options.q_lifetime || options.d_color)
    %start bin of second color
    red_start = 5000;

    blue = im_tcspc < red_start;
    red  = im_tcspc >= red_start;

    blue_tcspc      = im_tcspc(blue);
    blue_ISM_posx   = ISM_posx(blue);
    blue_ISM_posy   = ISM_posy(blue);

    red_tcspc      = im_tcspc(red);
    red_ISM_posx   = ISM_posx(red);
    red_ISM_posy   = ISM_posy(red);

    %crate blue ISM image
    [blue_ISM_img, blue_ISM_lin, ~]  = img_ps(blue_ISM_posx, blue_ISM_posy, ...
        s_pixl_x, s_pixl_y, ISM_binning);
    
    ISM_name_blue = append(ISM_name, '_blue');
    %plot and save blue
    img_plot(blue_ISM_img, hot, ISM_name_blue, 'blue ISM', 4, IM_R, ISM_binning, ...
        reso_line_ISM, 0, options.save_image);

    %crate red ISM image
    [red_ISM_img, red_ISM_lin, ~]  = img_ps(red_ISM_posx, red_ISM_posy, ...
        s_pixl_x, s_pixl_y, ISM_binning);
    
    ISM_name_red = append(ISM_name, '_red');
    %plot and save red
    img_plot(red_ISM_img, hot, ISM_name_red, 'red ISM', 4, IM_R, ISM_binning, ...
        reso_line_ISM, 0, options.save_image);

    if options.q_lifetime
        %lifetiem analysis
        [tail_t, tail_bin_l, tail_start_time, tcspc_t] = get_lifetime_bins(670, 5000, time_R, bin_factor);

        blue_pattern_tau = [1.5 3.6 inf];
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
        img_plot(medfilt2(lt_amp_img(:,:,1)), c_blue, lt_short_name{1}, shot_title, 4, IM_R, ISM_binning, reso_line_ISM, 0, 1);
        img_plot(medfilt2(lt_amp_img(:,:,2)), c_green, lt_long_name{1}, long_title, 4, IM_R, ISM_binning, reso_line_ISM, 0, 1);
    
        %lifetiem analysis red
        [tail_t, tail_bin_l, tail_start_time, tcspc_t] = get_lifetime_bins(5700, 9990, time_R, bin_factor);

        red_pattern_tau = [2.4 4.2 inf];
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
        img_plot(medfilt2(lt_amp_img(:,:,1)), c_yellow, lt_short_name{1}, shot_title, 4, IM_R, ISM_binning, reso_line_ISM, 0, 1);
        img_plot(medfilt2(lt_amp_img(:,:,2)), c_red, lt_long_name{1}, long_title, 4, IM_R, ISM_binning, reso_line_ISM, 0, 1);

    end
end

%% SOFI
if options.sofi
    ISM_SOFI_img = get_ISMSOFI(ISM_img, ISM_lin, im_time, head);

    img_plot(ISM_SOFI_img, spectrum, 'sofi ism', 'sofi ism', 1, IM_R, ISM_binning, reso_line_conf, 0, 1);

    SOFI_img = get_ISMSOFI(sum_img, sum_lin, im_time, head);

    img_plot(SOFI_img, spectrum, 'sofi sum', 'sofi sum', 1, IM_R, 1, reso_line_conf, 0, 1);
end

%% Additional figures for controle
if options.add_plt
    f = figure;
    ax = axes(f);
    histogram(ax, im_tcspc, 1:max_bin, 'EdgeAlpha',0)
    set(ax,'YScale', 'log')
    xlabel('tcspc bin')
    ylabel('count')
    
    if options.t_lifetime
        [count, ~]   = histcounts(im_tcspc, tcspc_t);
        
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
        
        amp1 = sum(lt_amp_img(:,:,1), 'all');
        amp2 = sum(lt_amp_img(:,:,2), 'all');
        amp3 = sum(lt_amp_img(:,:,3), 'all');
        offset = sum(lt_amp_img(:,:,4), 'all');
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
    end
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