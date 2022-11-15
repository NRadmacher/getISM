% Plot intensity distribution on spad detector for different z 
clear
close all
clc


names = dir('C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\210922\solution_OG_*.ptu');

lenght = size(names,1);

h = waitbar(0,'Lifetime fit');
for i = 1:lenght
    name = append(names(i).folder,'\', names(i).name);
    [im_chan,~,~,~,~, head] = read_ISM(name);
    
    title = head.ImgHdr_Z0;
    pixel_img   = mHist(im_chan,0:22);
    hex_plot(pixel_img,hot,title)
    
    waitbar(i/lenght,h)
end

close(h);