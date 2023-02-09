function [chan,tcspc,sync,head] = read_FCS(name)
%READ_ISM read 5d FLISM Data and header form PTU file 
%   Detailed explanation goes here
fname = name;

% head  = PTU_Read_Head(fname);

[sync, tcspc, chan, ~, ~, ~, head] = PTU_Read_old(name);
time_R = mean(head.MeasDesc_Resolution);
head.max_bin = round(1./head.TTResult_SyncRate /time_R);

end

