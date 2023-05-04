%% Clear shit  up
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F)
clear
% close all
clc

%% Parameters for analysis

options.lifetime    = 0;
options.s_lifetime  = 0;
options.d_lifetime  = 0;
options.d_color     = 0;
options.t_lifetime  = 0;
options.q_lifetime  = 0;

options.ISM_binning = 1.0;
options.WF_binning = 1.062;
options.ISM_sampling = 2;

% ymin xin ;ymax xmax
options.conf_rio    = [60 100; 120 160];
options.wf_rio      = [65 100; 125 160];
options.ISM_rio     = [120 200; 240 320];

options.reso_line_conf  = [[40 40]; [16 56]];
options.reso_line_wf    = [[46 46]; [16 56]];
options.reso_line_ISM   = [[74 74]; [30 110]];
options.reso_line_frw   = [[74 74]; [30 110]];

options.frw         = 1;
options.sofi        = 0;

options.add_plt     = 0;
options.save_image  = 1;
options.plot_reso   = 1;

%% PTU files with data, DC, IRF
fname   = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230424_ISM_STORM\gatta_quant_beadsR_50nm_pix_1.ptu';
dcname  = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230418_ISM_STORM\beads_calibration_PSF_200mm_laser_off6.ptu';
irfname = 'D:\PHD\Data\2022\220816\irf_ex470nm_005.ptu';

files = dir('C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230424_ISM_STORM\z_stack2\gatta_quant_beadsR_50nm_pix_*.ptu');

%% run
get_ISM(fname, dcname, irfname, options);

% for i = 1:size(files,1)
%     name = append(files(i).folder,'\', files(i).name);
%     get_ISM(name, dcname, irfname, options);
%     close all
% end
