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
options.d_color     = 1;
options.t_lifetime  = 0;
options.q_lifetime  = 0;

options.frw         = 0;
options.sofi        = 0;

options.add_plt     = 0;
options.save_image  = 1;
options.plot_reso   = 0;

%% PTU files with data, DC, IRF

fname   = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230413_CMs_fixed\CM_fixed_20kPa_ACTN2_citrine_Paxilin_SPAD_022.ptu';
dcname  = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\230314_CMs_life\dc_spad_001.ptu';
irfname = 'D:\PHD\Data\2022\220816\irf_ex470nm_005.ptu';



%% run

get_ISM(fname, dcname, irfname, options);
