function get_ISM(fname, dcname, irfname, caliname, options)

arguments
    fname char             
    dcname char           
    irfname char
    caliname char
    options (1,1) struct
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

%Image title and name from file name
tmp_name        = strsplit(fname, '\');
if contains(tmp_name{end-1},".")
    ws_name         = strsplit(tmp_name{end-1},'.');
    folder_name     = append(tmp_name{end-2},' ',ws_name{1});
else
    folder_name     = tmp_name{end-1};
end
date            = strsplit(folder_name, '_');
date            = date{1};
img_name        = tmp_name{end};
img_name        = strsplit(img_name, '.');
img_name        = img_name{end-1};
img_name        = append(date,' ',img_name);
ISM_name        = append(img_name, ' ISM');
ISM_fw_name     = append(img_name, ' ISM Fourier-reweighted');
PSF_name         = append(img_name, ' psf');
colf_name       = append(img_name, ' no ISM');
wf_name         = append(img_name, '  WF');
LT_name         = append(img_name, ' lt');

%read ISM data from .ptu file
% [im_chan,im_tcspc,im_posy,im_posx,im_time,head] = read_ISM(fname);

[im_time, im_tcspc, im_posx, im_posy, im_chan, head] = ScanRead(fname);

im_posx = double(im_posx);
im_posy = double(im_posy);

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

%shift vectors from file negativ sign is already included
calib = open(caliname);
sv = calib.shiftVector;
n_pixl = size(sv,2);


%ISM Binning UGO[1.02, 1.018] PQ[1.04] UGO220302[1.01] 20nm[1.0075] g4mix
%[1.01/1.0075]pmt
ISM_binning = options.ISM_binning;

%temporal resolution of TCSPC in seconds
time_R = mean(head.MeasDesc_Resolution);


if strcmp(head.CreatorSW_Name, 'SymPhoTime 64')
    %timer per pixel in total in seconds
    pix_time = head.ImgHdr_MaxFrames * head.ImgHdr_TimePerPixel/1e3;
else
    pix_time = head.ImgHdr_TimePerPixel/1e3;
end

max_bin = 1200;   
%max_bin = head.max_bin;

% TCSPC binning factor
bin_factor = 1;

% minimum numbers of photons to calculate lifetime
lt_cut_off = 10;

% tcspc tail_start in tcspc bin number !!FIX NEEDED!!
% tail_start = indMax + fwhm; diode [600/700,275], TiSa [250], SEPIA [1700]
tail_start = 130;

% tcspc tail_end in tcspc bin number
tail_end = max_bin - 20;

[tail_t, tail_bin_l, tail_start_time, tcspc_t] = get_lifetime_bins(tail_start, tail_end, time_R, bin_factor);

% Max Lifetime in ns
max_lt = 4;

% Normalised monoexponetial decay with background
pfun_monoexp = @(tau,x,dt) dt(:).*exp(-x(:)./tau); 
pfun_monoexpBG = @(tau,b,x,dt)b./numel(x(:))+(1-b).*pfun_monoexp(tau,x,dt)./sum(pfun_monoexp(tau,x,dt),1);

hex_plot(histcounts(im_chan,n_pixl), hot);
%% seperate frames
% 
% binning = 1;
% binFrame = 500;
% pix_time = 0;
% 
% sFrame = 1+(binFrame-1)*binning;
% eFrame = binFrame*binning;
% 
% 
% sTime = sFrame * head.ImgHdr_FrameTime;
% etime = eFrame * head.ImgHdr_FrameTime;
% ind = im_time<etime & im_time>sTime;%9.5;
% 
% im_tcspc        = im_tcspc(ind);
% im_chan         = im_chan(ind);
% im_posx         = im_posx(ind);
% im_posy         = im_posy(ind);
% im_time         = im_time(ind);

% s_pixl_y = head.ImgHdr_LineNum;
% s_pixl_x = head.ImgHdr_PixNum/head.ImgHdr_LineNum;

%% Confocal Image
%generate confocal image
[sum_img, sum_lin, ~] = img_ps(im_posx, im_posy, s_pixl_x, s_pixl_y,1);

r = figure('Visible','off');
r_ax = axes('Parent', r);

%plot and save
sum_img_dc = max(sum_img - sum(dc) * pix_time, 0);
img_plot( sum_img, hot, colf_name, ...
    'confocal', options.sb_lenght, IM_R, options.conf_rio, options.reso_line_conf, ...
    options.plot_reso, options.save_image, r, r_ax);

%% Wide Filed Image

if options.wf
    fprintf('Wide Filed reasigment ...');
    
    sv_file_name = 'SPAD_shift_vectors_wf.m';
    sv_file = matfile(sv_file_name);
    sv = (sv_file.sv_wf);
    shift_x     = sv(1, im_chan+1).';
    shift_y     = sv(2, im_chan+1).';
    
    %apply Wide Field reassigment vektor
    WF_posx    = im_posx + shift_x./(IM_R/0.05);
    WF_posy    = im_posy + shift_y./(IM_R/0.05);
    
    clear shift_y shift_x;
    fprintf('Done!\n');
    %generate confocal image
    [wf_img, ~, ~] = img_ps(WF_posx, WF_posy, s_pixl_x, s_pixl_y, options.WF_binning);
    
    WF_R = IM_R / options.WF_binning;
    
    %plot and save
    wf_img_dc = max(wf_img - sum(dc) * pix_time, 0);
    img_plot( wf_img_dc, hot, wf_name, ...
        'wide field', options.sb_lenght, WF_R, options.wf_rio, options.reso_line_wf, ...
        options.plot_reso, options.save_image, r, r_ax);
end

%% ISM reasigment
if options.ISM
    fprintf('ISM reasigment ... ');

    shift_x     = sv(1, im_chan+1).';
    shift_y     = sv(2, im_chan+1).';
    
%     [optBinning] = getISMbinning(s_pixl_y,sv, 5);
    
    ISM_binning = 1;
    % ISM_binning = options.ISM_binning;
    %apply ISM reassigment vektor
    ISM_posx    = im_posx + shift_x.*(0.05/IM_R);
    ISM_posy    = im_posy + shift_y.*(0.05/IM_R);
    
    clear shift_y shift_x;
    fprintf('Done!\n');
end
%% ISM Image
if options.ISM
    %crate ISM image
    [ISM_img, ISM_lin, ~]  = img_ps(ISM_posx, ISM_posy, s_pixl_x, s_pixl_y, ISM_binning);
    if options.ISM_sampling
        [ISM_img, ~,~] = FourierUpsampling(ISM_img, 2);
        ISM_R = IM_R / 2 / ISM_binning;
        ISM_pix_time = pix_time / 2;
    else
        ISM_R = IM_R / ISM_binning;
        ISM_pix_time = pix_time;
    end
    
    %calculate ISM darkcount: Dc is linear and every ISM pixel gets count from
    %23 pixels, either the same pixel in the sampel or a shifted on. But always
    %23. Thus darkcount = sum(dc)
    
    %plot and save
    img_plot(max(ISM_img - sum(dc) * ISM_pix_time, 0), hot, ISM_name, ...
        'ISM', options.sb_lenght, ISM_R, options.ISM_rio, options.reso_line_ISM, ...
        options.plot_reso, options.save_image, r, r_ax);
end
%% Fourier-reweighted ISM

if options.frw
    if options.ISM_sampling
        FW_R = ISM_R;
    else
        %double the sampling points to get reweighting at twice the k-vectors
        FW_R = ISM_R / 2;
    end
    
    %jörg PSF fit
%     [over, int] = NielsPSFFit(sum_img_dc);
    %calculate Fourier-reweighted ISM image
%     [W_ISM_img1, t_psf] = f_reweighting( max(ISM_img - sum(dc) * ISM_pix_time, 0), FW_R, calib.over, calib.psf);
    W_ISM_img1 = ISM_frw(ISM_img,calib.psf);
    %plot and save fr ISM
    img_plot(W_ISM_img1, hot, ISM_fw_name, 'Fourier reweighted ISM', ...
        options.sb_lenght, ISM_R, options.ISM_rio, options.reso_line_frw, ...
        options.plot_reso, options.save_image, r, r_ax);

    if options.add_plt
        %iamge size 
        t_psf = calib.psf;
        nx = -options.conf_rio(1,1)+options.conf_rio(2,1);
        ny = -options.conf_rio(1,2)+options.conf_rio(2,2);
        [Nx,Ny] = size(t_psf);
        if Nx<nx || Ny<ny
            nx = Nx;
            ny = Ny;
        end
        t_psf = t_psf(floor((Nx-nx)/2)+(1:nx), floor((Ny-ny)/2)+(1:ny));

        img_plot(t_psf, hot, PSF_name, 'confocal PSF', ...
        0.5, IM_R, 0, [[17 17]; [16-10 16+10]], ...
        options.plot_reso, options.save_image, r, r_ax);
    end
end
%% Plot reso

if options.plot_reso
    set(r, 'Visible', 'on')
    grid(r_ax, 'on');
    box(r_ax, 'on');
%     xlim(r_ax,[0.5 2]);
    lgd = legend(r_ax, 'confocal', '','ISM','','FRW ISM');
    lgd.Location = 'best';
    file_name = append(img_name,'_reso.pdf');
    exportgraphics(r_ax, file_name,'Resolution',600)
end
%% Single Expoential lifetime
if options.s_lifetime
    %get pixel lifetimes via MLE pattern matching
    lt_img = get_single_lifetime(ISM_img,ISM_lin,...
        im_tcspc = im_tcspc, tail_t = tail_t, tail_bin_l = tail_bin_l,...
        tail_start_time = tail_start_time, tcspc_t = tcspc_t,...
        max_lt = max_lt, lt_cut_off = lt_cut_off, fname = LT_name, ...
        ana_plt = options.add_plt);

    %plot pxel wise lt values scaled with image intensity
    lt_img_plot(lt_img, ISM_img, c_map, lt_cut_off, options.lt_range, ...
        'Lifetime', LT_name, 4, ISM_R, ISM_binning, options.save_image)   
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
    
    % tail time that starts at 0
    xx = tail_t-tail_start_time;

    % tcspc histogram
    [count, ~]   = histcounts(im_tcspc, tcspc_t);
    [~,tau] = double_exp(xx.',count.');
    pattern_tau = tau;
    pattern = pfun_monoexpBG(pattern_tau,0,tail_t-tail_start_time,tail_bin_l);

    % Struckture and Fluorophore name
    names = {'tubulin', 'actin'; 'AF 647', 'JF 646'};

    lt_short_name   = compose('%s lt short PIRL %0.1f ns', ISM_name, pattern_tau(1)); %1.37
    lt_long_name    = compose('%s lt long PIRL %0.1f ns', ISM_name, pattern_tau(2)); % 3.0

    shot_title = append(names{1,1}, ': ', names{2,1});
    long_title = append(names{1,2}, ': ', names{2,2});
    
    %unmix the two patterns
    [lt_amp_img] = get_bi_lifetime(ISM_img,ISM_lin,im_tcspc = im_tcspc, ...
        tail_t = tail_t, tail_bin = tail_bin_l, tail_start_time = tail_start_time, ...
        tcspc_t = tcspc_t, pattern = pattern, pattern_tau = pattern_tau, ...
        name = names, fname = ISM_name, ana_plt = options.add_plt);
    
    %plot individual strucktures
    img_plot(lt_amp_img(:,:,1), c_red, lt_short_name{1}, shot_title, 4, ...
        ISM_R, options.ISM_rio, options.reso_line_ISM, 0, 1, r, r_ax);
    img_plot(lt_amp_img(:,:,2), c_green, lt_long_name{1}, long_title, 4, ...
        ISM_R, options.ISM_rio, options.reso_line_ISM, 0, 1, r, r_ax);

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

    img_plot(ISM_SOFI_img, spectrum, 'sofi ism', 'sofi ism', 1, IM_R, ISM_binning, options.reso_line_ISM, 0, 1);

    SOFI_img = get_ISMSOFI(sum_img, sum_lin, im_time, head);

    img_plot(SOFI_img, spectrum, 'sofi sum', 'sofi sum', 1, IM_R, 1, options.reso_line_ISM, 0, 1);
end
%% Additional figures for controle
if options.add_plt
    f = figure;
    ax = axes(f);
    histogram(ax, im_tcspc, 1:max_bin, 'EdgeAlpha',0)
    set(ax,'YScale', 'log')
    xlabel('tcspc bin')
    ylabel('count')

    plot_pixeldecay(im_tcspc,im_chan,32)
    
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
end