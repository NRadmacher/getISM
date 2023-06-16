%% Clear shit  up
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F)
clear
close all
clc

%% Parameters for analysis

options.lifetime    = 0;
options.s_lifetime  = 1;
options.d_lifetime  = 0;
options.d_color     = 0;
options.t_lifetime  = 0;
options.q_lifetime  = 0;

%1.02
options.ISM_binning = 1.08000;
options.WF_binning = 1.062;
options.ISM_sampling = 0;

% ymin xin ;ymax xmax
options.conf_rio    =0;%[75 108; 115 148];%[60 100; 120 160];
options.wf_rio      =0;%[65 100; 125 160];
options.ISM_rio     =0;%[150 210; 230 290];%[120 200; 240 320];

% options.ISM_rio     = [60 100; 120 160];

% options.reso_line_conf  = [[40 40]; [16 56]];
% options.reso_line_ISM   = [[74 74]; [30 110]];
% options.reso_line_frw   = [[74 74]; [30 110]];

options.reso_line_conf  = [[32 32]; [11 31]];
options.reso_line_ISM   = [[64 64]; [20 60]];
options.reso_line_frw   = [[64 64]; [20 60]];

% options.reso_line_ISM   = [[36 36]; [15 55]];
% options.reso_line_frw   = [[36 36]; [15 55]];

options.lt_range = [1.5 1.8];

options.wf          = 0;
options.frw         = 0;
options.sofi        = 0;

options.add_plt     = 0;
options.save_image  = 1;
options.plot_reso   = 0;

%% PTU files with data, DC, IRF
fname   = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230615_ISM_STORM\neurons_4pfa_psd95_syt1_alexa647_atto655_PBS_power6p5_OD3_roi1_3.ptu';
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
