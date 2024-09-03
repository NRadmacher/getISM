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
ind(12) = 0;
ind(5) = 0;
ind = logical(ind);

%% 

calib230711 = open("ismCallibration230711_ISM_STORM.mat");
calib240807 = open("ismCallibration240807_ism.mat");
calib240823 = open("ismCallibration240823_ismAligment.mat");
calib240828 = open("ismCallibration240828_ismAligment.mat");
calib240830 = open("ismCallibration240830_ism.mat");
calib240902 = open("ismCallibration240902_ismAligment.mat");

norm230711 = vecnorm(calib230711.shiftVector,2,1).*calib230711.pixSize;
norm240807 = vecnorm(calib240807.shiftVector,2,1).*calib240807.pixSize;
norm240823 = vecnorm(calib240823.shiftVector,2,1).*calib240823.pixSize;
norm240828 = vecnorm(calib240828.shiftVector,2,1).*calib240828.pixSize;
norm240830 = vecnorm(calib240830.shiftVector,2,1).*calib240830.pixSize;
norm240902 = vecnorm(calib240902.shiftVector,2,1).*calib240902.pixSize;

alpha = 1+708/640;
alpha = 2;
M = mean(dist(ind)./(alpha*norm230711(ind)),"omitnan");
disp(M)
M = mean(dist(ind)./(alpha*norm240807(ind)),"omitnan");
disp(M)
M = mean(dist(ind)./(alpha*norm240823(ind)),"omitnan");
disp(M)
M = mean(dist(ind)./(alpha*norm240828(ind)),"omitnan");
disp(M)
M = mean(dist(ind)./(alpha*norm240830(ind)),"omitnan");
disp(M)
M = mean(dist(ind)./(alpha*norm240902(ind)),"omitnan");
disp(M)



