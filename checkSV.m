%% Clear shit  up
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F)
clear
close all
clc


%% spad dinesions
a = 1;
b = 2;
c = sqrt( (1*sqrt(3)/2)^2 + (1/2)^2 );
d = sqrt( (1*sqrt(3)/2)^2 + (3/2)^2 );
e = sqrt( (2*sqrt(3)/2)^2 + (0)^2 );
f = sqrt( (2*sqrt(3)/2)^2 + (1)^2 );
g = sqrt( (2*sqrt(3)/2)^2 + (2)^2 );
dist = [g,f,e,f,g,d,c,c,d,b,a,0,a,b,d,c,c,d,g,f,e,f,g].*23;
ind = ones(1,23);
ind(12) = 1;
ind(5) = 1;
ind = logical(ind);

%% 

calib = open("ismCallibration240927 lifetime GroupMeas_2 gattaBeadsRed_18_z0,80µm_1.mat");

norm = vecnorm(calib.shiftVector,2,1).*calib.pixSize;

alpha = 1+(708/640);
% alpha = 2;

M = mean(dist(ind)./(alpha*norm(ind)),"omitnan");
disp(M)


