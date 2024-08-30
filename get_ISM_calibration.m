%claculate ISM shift Vercores and PSF

% fname = 'W:\Niels\Messungen_Daten\240426_THG_ISM\GaAs_edge_THG_008.ptu';
% dcname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230711_ISM_STORM\gattaQ_beads_atto647n_1p9mW_OD2_50nm-pix_dark5.ptu';

% fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230711_ISM_STORM\gattaQ_beads_atto647n_1p9mW_OD2_50nm-pix_2.ptu';
% dcname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230711_ISM_STORM\gattaQ_beads_atto647n_1p9mW_OD2_50nm-pix_dark5.ptu';

fname = 'W:\Niels\Messungen_Daten\240830_ism\lifetime.sptw\gattaBeadsRed_3.ptu';
dcname = 'W:\Niels\Messungen_Daten\240830_ism\lifetime.sptw\dc_2.ptu';

%% Microscope parameters for psf fit
%numerical aperture
NA = 1.49;
%focal distance of objective [µm]
fd = 1800;
%excitation wavelenght [µm]
lamex = 0.640;

%% Other parameters
%total maginfication at array
M = 150;

%% Load data
tmp_name        = strsplit(fname, '\');
if contains(tmp_name{end-1},".")
    folder_name     = tmp_name{end-2};
else
    folder_name     = tmp_name{end-1};
end

%read ISM data from .ptu file
[~, ~, im_posx, im_posy, ~, head] = ScanRead(fname);

im_posx = double(im_posx);
im_posy = double(im_posy);

ind             = im_posy<=head.ImgHdr_PixY;
im_posx         = im_posx(ind);
im_posy         = im_posy(ind);
im_posx         = im_posx * head.ImgHdr_PixX;

%number of pixels in recorded image
s_pixl_x    = head.ImgHdr_PixX;
s_pixl_y    = head.ImgHdr_PixY;

%Dark count in count per second
[dc, ~] = get_DC(dcname,0);

%time per pixel
if strcmp(head.CreatorSW_Name, 'SymPhoTime 64')
    %timer per pixel in total in seconds
    pix_time = head.ImgHdr_MaxFrames * head.ImgHdr_TimePerPixel/1e3;
else
    pix_time = head.ImgHdr_TimePerPixel/1e3;
end

calibration = struct();

calibration.dc = dc;
%% calculate confocal psf

[sum_img, ~, ~] = img_ps(im_posx, im_posy, s_pixl_x, s_pixl_y,1);
sum_img_dc = max(sum_img - sum(dc) * pix_time, 0);

[over, int, im, xx, yy] = ismPSFFit(sum_img,...
    head.ImgHdr_PixResol,NA,fd,lamex);

calibration.psf     = int;
calibration.over    = over;
calibration.NA      = NA;
calibration.fd      = fd;
calibration.lamex   = lamex;
calibration.PSFfunc = @PSF;
%% generate shift Vectors

[shiftVector, save_name] = getShiftVectors(fname, dcname);

calibration.shiftVector = shiftVector;
calibration.M           = M;
calibration.svName      = save_name;
calibration.fileName    = fname;
calibration.dc          = dc;
calibration.pixSize     = head.ImgHdr_PixResol;

saveTo = append('ismCallibration', folder_name, '.mat');

save(saveTo,"-struct", "calibration",'-v7.3')