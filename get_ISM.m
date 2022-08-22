%% Clear shit  up
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F)
clear
close all
clc

lifetime    = 0;
deconv      = 0;
sofi        = 1;
%% Loadind data
%rgb values for color map black,blue,cyan,green,yellow,orange?,red,magenta
sp1     = 1:255/7:256;
sp2     = [[0 0 0]; [0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]; [1 0 1]];
br1     = 1:255/5:256;
br2     = [[0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]];
grl     = 1:255/3:256;
gr2     = [[0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]];
gl      = 1:255/1:256;
g2      = [[0 0 0]; [0.454325 0.806075 0.332481]];
% g2      = [[0 0 0];[0 0.9686 0]];
rl      = 1:255/1:256;
r2      = [[0 0 0]; [1.000000 0.630093 0.501411]];
% r2      = [[0 0 0]; [0.9961 0 0]];
bl      = 1:255/1:256;
b2      = [[0 0 0]; [0.000000 0.783196 1.000000]];
% b2      = [[0 0 0]; [0.9608 1 0]];

lambda  = 1:256;

spectrum    = interp1(sp1, sp2, lambda);
c_greenred  = interp1(grl, gr2, lambda);
c_bluered   = interp1(br1, br2, lambda);
c_green     = interp1(gl, g2, lambda);
c_red       = interp1(rl, r2, lambda);
c_blue      = interp1(bl, b2, lambda);
c_map       = cmap_isoluminant75;

fname   = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220726\qdots_em605nm_012.ptu';
dcname  = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\220106\cd_001.ptu';
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
[im_chan,im_tcspc,im_posy,im_posx,im_time, head] = read_ISM(fname);

ind             = im_posy<=head.ImgHdr_PixY;
im_chan         = im_chan(ind);
im_posx         = im_posx(ind);
im_posy         = im_posy(ind);
im_posx         = im_posx * head.ImgHdr_PixX;
im_time         = im_time./head.TTResult_SyncRate; % photon arrival in seconds

%remove time delay due to unsyncronised Multi Harps
% im_tcspc = remove_MHH_offset(im_tcspc,im_chan, head.max_bin, 100);

%Dark count in count per second
[dc, bin_dc] = get_DC(dcname,0);

%get IRF to set tail for tailfit
% [fwhm, indMax, ~, ~,irf_tcspc] = get_IRF(irfname, 1);
%% Parameters and magic numbers(please fix) AND FIX NAMEN FÜR TCSPC

%number of pixels in recorded image
s_pixl_x    = head.ImgHdr_PixX;
s_pixl_y    = head.ImgHdr_PixY;

%scan resolution step per pixel in µm
IM_R = head.ImgHdr_PixResol;

%ISM Binning UGO[1.02] PQ[1.04] UGO220302[1.01] 20nm[1.0075]
ISM_binning = 1.01;

%temporal resolution of TCSPC in seconds
time_R = mean(head.MeasDesc_Resolution);

%number of events
n_events = numel(im_posx);

%max number of tcspc bins %%CARE
max_bin = head.max_bin;

% TCSPC binning factor
bin_factor = 20;

% minimum numbers of photons to calculate lifetime
lt_cut_off = 0;

% tcspc tail_start in tcspc bin number !!FIX NEEDED!!
% tail_start = indMax + fwhm;
tail_start = 760;

% tcspc tail_end in tcspc bin number
tail_end = max_bin - 100;

% tcspc tail lenght in tcspc bin number
tail_l = tail_end - tail_start;

% TCSPC bin length in ns
tcspc_bin_l = time_R * 1e9;

% combin binning for lifetime fitting
% Lifetime bin lenght in ns
lt_bin_l = tcspc_bin_l * bin_factor;

% Lifetime tail 
lt_start    = floor(tail_start / bin_factor);
lt_end      = floor(tail_end / bin_factor); %floor or ciel?

% TCSPC tail lenght in ns
% tcspc_tail_l = (tail_end - tail_start) * tcspc_bin_l;
tcspc_tail_l = (lt_end - lt_start) * lt_bin_l;

% TCSPC tail start in ns
tcspc_start = tail_start * tcspc_bin_l;

% Max Lifetime in ns
max_lt = 8;

% Normalised monoexponetial decay with background
pfun_monoexp = @(tau,x,dt) dt(:).*exp(-x(:)./tau); 
pfun_monoexpBG = @(tau,b,x,dt)b./numel(x(:))+(1-b).*pfun_monoexp(tau,x,dt)./sum(pfun_monoexp(tau,x,dt),1);

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
    tail_bin    = tcspc_bin_l*tcspc_bin(1:end-1);
    tail_t      = tcspc_bin_l*tcspc_t(1:end-1);

else
    tcspc_bin   = bin_factor;
    tcspc_t     = tail_start:bin_factor:tail_end;

    % and in ns for tail fit
    tail_bin    = tcspc_bin_l*tcspc_bin;
    tail_t      = tcspc_bin_l*tcspc_t(1:end-1);
end

%% Dark Count

bin_dc = sum(bin_dc, 1);

bin_dc = movsum(bin_dc, [0 bin_factor-1]);

bin_dc = bin_dc(tcspc_t);

%% ISM reasigment
fprintf('ISM reasigment\n');
%shift vectors from file negativ sign is already included
sv_file = matfile('shift_vector.m');
sv = sv_file.shift_vector;
shift_x     = sv(1, im_chan+1).';
shift_y     = sv(2, im_chan+1).';

%apply ISM reassigment vektor
ISM_posx    = im_posx - shift_x;
ISM_posy    = im_posy - shift_y;

%% Confocal Image
%generate confocal image
[sum_img, sum_lin, ~] = img_ps(im_posx, im_posy, s_pixl_x, s_pixl_y,1);

%remoce dark count
sum_img = max(sum_img - sum(dc) * head.ImgHdr_PixelTime, 0);

%plot and save
reso_line_conf = [[105 105]; [30 65]];
img_plot(sum_img, spectrum, colf_name, 'confocal', 4, IM_R, reso_line_conf, 0, 1);

%% ISM Image
%crate ISM image
[ISM_img, ISM_lin, ISM_size]  = img_ps(ISM_posx, ISM_posy, s_pixl_x, s_pixl_y, ISM_binning);

%calculate ISM darkcount: Dc is linear and every ISM pixel gets count from
%23 pixels, either the same pixel in the sampel or a shifted on. But always
%23. Thus darkcount = sum(dc)
% ISM_img  = max(ISM_img - sum(dc) * head.ImgHdr_PixelTime, 0);

%plot and save
reso_line_ISM = [[107 107]; [30 65]];
img_plot(max(ISM_img - sum(dc) * head.ImgHdr_PixelTime, 0), spectrum, ISM_name, 'ISM', 4, IM_R, reso_line_ISM, 0, 1);

%% manuel lucy Richerson decon
if(deconv)
    PSF_file = matfile('PSF.m');
    psf      = PSF_file.im;
    psf      = psf/sum(psf, 'all');
    
    EID_file = matfile('EID.m');
    eid      = EID_file.im;
    eid      = eid/sum(eid, 'all');
    
    manuel_decon(ISM_img, eid, ...
                c_map = spectrum, s_name = ISM_docn_name, ...
                t_name = 'ISM + deconvolution', ...
                sb_lenght = 1,IM_R = IM_R,...
                reso = 0, save = 1);
end

%% Lifetime Image
if (lifetime)

    interval        = ceil(1.01* max(max(ISM_img)));
    n_pixel_ISM     = prod(ISM_size);
    ISM_lt          = zeros(n_pixel_ISM, numel(tail_t));
    ISM_lt_amp      = zeros(n_pixel_ISM, 4);
    ISM_lt_rgb      = zeros([ISM_size 3]);

    % generate decay patterns from FL selecion[1.37 2.36 3.0]
    pattern_tau = [0.56 2.3 3.32 inf];
    pattern = pfun_monoexpBG(pattern_tau, 0, tail_t-tcspc_start, tail_bin);

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
%                max(count - 0.8*bin_dc(1:end-1).*head.ImgHdr_PixelTime, 0);
%                count        = count./sum(count);
%                ISM_lt(i,:)  = count;
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
% test tail fit
%     
%     ISM_img = ones(ISM_size)*20;
%     lin_tau = ceil(linspace(0.01, max_lt, n_pixel_ISM));
% 
%     ISM_lt = pfun_monoexpBG(lin_tau, 0.01, tail_t-tcspc_start, tail_bin);
% 
% %     noise = wgn( numel(tail_t), n_pixel_ISM, 0);
%     noise = normrnd(0, 1e-3,numel(tail_t), n_pixel_ISM);
%     ISM_lt = ISM_lt + noise;
% 
%     figure
%     hold on
%     plot( ISM_lt(:,50000) )
%     plot( ISM_lt(:,10) )
%     plot( ISM_lt(:,100000) )

    %tail fit via pattern matching
%     [ISM_lt_img,~]  = lt_patternMatching(ISM_lt, tail_t-tcspc_start, tail_bin, max_lt);
% 
% %     err = ISM_lt_img.' - lin_tau;
% % 
% %     figure
% %     swarmchart(lin_tau, err, 'XJitterWidth', 0.3)
% 
%     %plot and save image
%     ISM_lt_img  = reshape(ISM_lt_img, ISM_size);
%     lt_img_plot(ISM_lt_img, ISM_img, c_map, lt_cut_off, [0.1 max_lt], 'Lifetime', LT_name, 2, IM_R, 1)
% 
%     h = figure;
%     ax = axes(h);
%     histogram(ISM_lt_img(ISM_lt_img > 0.1 & ISM_lt_img < max_lt),linspace(0.01,max_lt,100),'Normalization','count')
% 
%     file_name = append(LT_name, '_dist', '.png');
%     exportgraphics(ax, file_name,'Resolution',600)

    ISM_lt_amp_img = reshape(ISM_lt_amp, [ISM_size(1), ISM_size(2), 4]);
    img_plot(ISM_lt_amp_img(:,:,1), c_green, ISM_lt_short_name{1}, 'Cy2: SYT 1', 4, IM_R, reso_line_ISM, 0, 1);
    img_plot(ISM_lt_amp_img(:,:,2), c_red, ISM_lt_middle_name{1}, 'OG: GFAP', 4, IM_R, reso_line_ISM, 0, 1);
    img_plot(ISM_lt_amp_img(:,:,3), c_blue, ISM_lt_long_name{1}, 'Alexa PSD95', 4, IM_R, reso_line_ISM, 0, 1);

    ISM_lt_rgb(:,:,1) = ISM_lt_amp_img(:,:,1)./max(ISM_lt_amp_img(:,:,1),[],'all');%red
    ISM_lt_rgb(:,:,2) = ISM_lt_amp_img(:,:,3)./max(ISM_lt_amp_img(:,:,3),[],'all');%green
    ISM_lt_rgb(:,:,3) = ISM_lt_amp_img(:,:,2)./max(ISM_lt_amp_img(:,:,2),[],'all');%blue

    comb_name = compose('%s multicolor %0.1f %0.1f %.01f', LT_name, pattern_tau(1), pattern_tau(2), pattern_tau(3));

    rgb_img_plot(ISM_lt_rgb, comb_name{1}, 'Red: SYT 1, Green: PSD95, Blue: GFAP', 4, IM_R, 1)

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

%% SOFI
if (sofi)

    [SOFI_img, ISM_SOFI_img] = get_SOFI(sum_img,sum_lin, im_time, im_chan, sv, head);

    img_plot(ISM_SOFI_img, spectrum, 'sofi ism', 'sofi ism', 1, IM_R, reso_line_conf, 0, 0);
    img_plot(SOFI_img, spectrum, 'sofi', 'sofi', 1, IM_R, reso_line_conf, 0, 0);
end

%% Additional figures for controle

figure
histogram(im_tcspc, 1:max_bin)

[count, edges]   = histcounts(im_tcspc, tcspc_t);

[tripple_amp,tau] = tripple_exp(tail_t.'-tcspc_start,count.');
disp(tau)
disp(tripple_amp)

xx = tail_t-tcspc_start;
count = count./sum(count);
[over_all_lt,over_all_back]  = lt_patternMatching(count, tail_t-tcspc_start, tail_bin, max_lt);
disp(over_all_lt)

figure
plot(tail_t-tcspc_start, count, 'bo')
hold on

yy = pfun_monoexpBG(over_all_lt, over_all_back, tail_t, tail_bin);
plot(tail_t-tcspc_start, yy.','r-')

amp1 = tripple_amp(1);
amp2 = tripple_amp(2);
amp3 = tripple_amp(3);
offset = tripple_amp(4);
y4 =  amp1*exp(-xx/tau(1)) + amp2*exp(-xx/tau(2)) + amp3*exp(-xx/tau(3)) + offset;
y4 = y4./sum(y4);
plot(xx,y4,'g-')

amp1 = sum(ISM_lt_amp_img(:,:,1), 'all');
amp2 = sum(ISM_lt_amp_img(:,:,2), 'all');
amp3 = sum(ISM_lt_amp_img(:,:,3), 'all');
offset = sum(ISM_lt_amp_img(:,:,4), 'all');
y3 =  amp1*exp(-xx/pattern_tau(1)) + amp2*exp(-xx/pattern_tau(2)) + amp3*exp(-xx/pattern_tau(3)) + offset;
y3 = y3./sum(y3);
plot(xx,y3,'m-')

legend('data','pattern Matching', 'true tripple', 'good looking')
set(gca, 'YScale', 'log')
ylim([0.5*min(yy) inf])
xlabel('time [ns]')
ylabel('normalized count')



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

