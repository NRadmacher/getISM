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
confx = 52;
confy = 162;

ismx = 52;
ismy = 103;
options.reso_line_conf  = [[confx-20 confx+20]; [confy confy]];
options.reso_line_ISM   = [[ismx ismx]; [ismy-20 ismy+20]];
options.reso_line_frw   = options.reso_line_ISM;%[[323 323]; [216 256]];

options.lt_range = [0.1 5];

options.ISM         = 0;
options.wf          = 0;
options.frw         = 0;
options.sofi        = 0;

options.eps = 0.08;

options.add_plt     = 0;
options.save_image  = 1;
options.plot_reso   = 1;

options.pixShift = 0;

%% PTU files with data, DC, IRF
fname       = 'W:\Niels\Messungen_Daten\240927_ism\lifetime.sptw\GroupMeas_2\gattaBeadsRed_18_z0,80µm_1.ptu';
dcname      = 'W:\Niels\Messungen_Daten\240806_ism\lifetime.sptw\dc_WL_2.ptu';
irfname     = 'D:\PHD\Data\2022\220816\irf_ex470nm_005.ptu';
caliname    = 'ismCallibration240927 lifetime GroupMeas_2 gattaBeadsRed_18_z0,80µm_1.mat';

files       = dir('W:\Niels\Messungen_Daten\240927_ism\lifetime.sptw\GroupMeas_2\gattaBeadsRed*.ptu');

%% run
get_ISM(fname, dcname, irfname, caliname, options);

% for i = 1:size(files,1)
%     name = append(files(i).folder,'\', files(i).name);
%     try 
%         get_ISM(name, dcname, irfname,caliname, options);
%     catch error
%         disp(error.message)
%     end
%     close all
% end
