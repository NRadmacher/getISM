%%
%calculate theoretical sv and psf

%% Microscope parameters for psf fit
%numerical aperture
NA = 1.35;
%focal distance of objective [µm] for Olympus = 180(mm)/M
fd = 180/40*1000;
%excitation wavelenght [µm]
lamex = 1.040;

over = 1800;

%adjustment for shiftet back and forward strock
pixShift = 0;

% t = datetime('today','ConvertFrom','yyyymmdd');

%% Other parameters
%total maginfication at array
M = 270;

%scan pixel size [µm]
pix = 0.05;

calibration = struct();

dc = 1;

%% PSF
calibration.over    = over;
calibration.NA      = NA;
calibration.fd      = fd;
calibration.lamex   = lamex;
calibration.PSFfunc = @PSF;
calibration.psf     = PSF(NA,fd,lamex,pix,over);
disp("PSF fit done!")

%% generate shift Vectors

sv_wf = get_th_shift_vectors();

calibration.shiftVector = sv_wf;
calibration.M           = M;
calibration.svName      = '_thSV_SPAD23_SHG';
calibration.fileName    = 'noName';
calibration.dc          = dc;
calibration.pixSize     = pix;
disp("shiftvectors calculated!")


saveTo = append('ismCallibration', 'thForSHG', '.mat');
save(saveTo,"-struct", "calibration",'-v7.3')



