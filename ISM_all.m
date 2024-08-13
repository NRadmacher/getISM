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

%1.02
options.ISM_binning = 1.0125;%1.08;
options.WF_binning = 1.062;
options.ISM_sampling = 0;

% ymin xin ;ymax xmax
options.conf_rio    = [95 160; 130 195];%[60 100; 120 160];
options.wf_rio      = [65 100; 125 160];
options.ISM_rio     = [95 160; 130 195];%[120 200; 240 320];

options.conf_rio    = 0;
options.wf_rio      = 0;
options.ISM_rio     = 0;

options.sb_lenght   = 5;

options.reso_line_conf  = [[118 118]; [144 184]];
options.reso_line_ISM   = [[120 120]; [144 184]];
options.reso_line_frw   = [[120 120]; [144 184]];

options.lt_range = [0.1 5];

options.ISM         = 1;
options.wf          = 0;
options.frw         = 1;
options.sofi        = 0;

options.add_plt     = 0;
options.save_image  = 1;
options.plot_reso   = 0;

%% PTU files with data, DC, IRF
fname       = 'W:\Niels\Messungen_Daten\240801_ism\exAhi.sptw\otoferlin_Star635p_gold_ODp8_1.ptu';
dcname      = 'W:\Niels\Messungen_Daten\240806_ism\lifetime.sptw\dc_WL_2.ptu';
irfname     = 'D:\PHD\Data\2022\220816\irf_ex470nm_005.ptu';
caliname    = 'ismCallibration240807_ism.mat';

files       = dir('W:\Niels\Messungen_Daten\240801_ism\exAhi.sptw\*.ptu');

%% run
get_ISM(fname, dcname, irfname, caliname, options);

% for i = 1:size(files,1)
%     name = append(files(i).folder,'\', files(i).name);
%     get_ISM(name, dcname, irfname,caliname, options);
%     close all
% end
