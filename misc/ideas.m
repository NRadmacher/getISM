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


bin_dc = sum(bin_dc, 1);

bin_dc = movsum(bin_dc, [0 bin_factor-1]);

bin_dc = bin_dc(tcspc_t);

%% manuel lucy Richerson decon
if(deconv)
    PSF_file = matfile('PSF.m');
    psf      = PSF_file.im;
    EID_file = matfile('EID.m');
    eid      = EID_file.im;

    [eid, psf] = get_psf([0 0.8], 0.05);
    psf      = psf/sum(psf, 'all');
    eid      = eid/sum(eid, 'all');
    
    manuel_decon(max(ISM_img - sum(dc) * head.ImgHdr_PixelTime, 0), eid, ...
                c_map = hot, s_name = ISM_docn_name, ...
                t_name = 'ISM + deconvolution', ...
                sb_lenght = 1, IM_R = IM_R, pix_bin = ISM_binning, ...
                reso = 1, save = 1);
end

%% Lifetime Image

%TODO see if this makes sense
    test_img = ISM_lt_amp_img(:,:,1:3);
    test_img = ISM_lt_amp_img./repmat(max(ISM_lt_amp_img,[],[1 2]), [size(ISM_lt_amp_img,[1 2]) 1]);
    test_img = test_img(:,:,1:3);
    test_color =[g2(2,:); r2(2,:); b2(2,:)];
    test_color =[[1 0 1]; [0 1 1]; [0.96 1 0]];

    test_rgb = tensorprod(test_img,test_color,3,1);
    test_rgb = test_rgb./repmat(max(test_rgb,[],[1 2]), [size(test_rgb,[1 2]) 1]);
    rgb_img_plot(test_rgb, LT_name, 'Red: GFAP, Green: SYT, Blue: PSD95', 4, IM_R, 0)


%%
    % for sparse sample (at best SM) sqaring is the same a FW
    img_plot((max(ISM_img - sum(dc) * ISM_pix_time, 0).^2), hot, ISM_sq_name, ...
    'ISM real Raum quadrat', 1, ISM_R, ISM_binning, options.ISM_rio, options.reso_line_ISM, ...
    options.plot_reso, options.save_image);

%% JE lt fit
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
%% focus ISM

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