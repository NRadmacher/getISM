%% Clear shit  up
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F)
clear
close all
clc

%% Parameters for analysis

options.lifetime    = 0;
options.s_lifetime  = 0;
options.d_lifetime  = 0;
options.d_color     = 0;
options.t_lifetime  = 0;
options.q_lifetime  = 0;

options.ISM_binning = 1.02;

options.frw         = 0;
options.sofi        = 0;

options.add_plt     = 0;
options.save_image  = 1;
options.plot_reso   = 0;

%% PTU files with data, DC, IRF
fname   = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230421_ISM_STORM\gatta_bead_atto647N_50nm_pix_3.ptu';
dcname  = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230418_ISM_STORM\beads_calibration_PSF_200mm_laser_off6.ptu';
irfname = 'D:\PHD\Data\2022\220816\irf_ex470nm_005.ptu';

files = dir('C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230413_CMs_fixed\*.ptu');

%% run
get_ISM(fname, dcname, irfname, options);

% for i = 1:size(files,1)
%     name = append(files(i).folder,'\', files(i).name);
%     get_ISM(name, dcname, irfname, options);
%     close all
% end
