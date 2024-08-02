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
options.ISM_binning = 1;%1.08;
options.WF_binning = 1.062;
options.ISM_sampling = 0;

% ymin xin ;ymax xmax
options.conf_rio    = [95 160; 130 195];%[60 100; 120 160];
options.wf_rio      = [65 100; 125 160];
options.ISM_rio     = [95 160; 130 195];%[120 200; 240 320];

options.conf_rio    = 0;
options.wf_rio      = 0;
options.ISM_rio     = 0;

options.sb_lenght   = 1;

options.reso_line_conf  = [[84 84]; [24 74]];
options.reso_line_ISM   = [[82 82]; [24 74]];
options.reso_line_frw   = [[82 82]; [24 74]];

options.lt_range = [1.5 1.8];

options.ISM         = 0;
options.wf          = 0;
options.frw         = 0;
options.sofi        = 0;

options.add_plt     = 0;
options.save_image  = 1;
options.plot_reso   = 0;

%% PTU files with data, DC, IRF
fname = 'W:\Niels\Messungen_Daten\240612_ISM\cali.sptw\gattaBeadsR_49.ptu';
dcname = 'W:\Niels\Messungen_Daten\240612_ISM\cali.sptw\dcFastFrameRate_1.ptu';
irfname = 'D:\PHD\Data\2022\220816\irf_ex470nm_005.ptu';

files = dir('W:\Niels\Messungen_Daten\240612_ISM\cali.sptw\gattaBeadsR_*.ptu');

%% run
get_ISM(fname, dcname, irfname, options);

% for i = 1:size(files,1)
%     name = append(files(i).folder,'\', files(i).name);
%     get_ISM(name, dcname, irfname, options);
%     close all
% end
