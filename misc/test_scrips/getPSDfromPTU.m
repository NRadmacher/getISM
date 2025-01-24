


files = dir('W:\Niels\Messungen_Daten\250124_ism\perovskiteTransport.sptw\short_fastLT_1\TCMH_SNO2_no1_pulsed_10MHz_OD3*.ptu');

f = figure;
ax = axes(f);

for i = 1:size(files,1)
    name = append(files(i).folder,'\', files(i).name);
    disp([i size(files,1)])
    try 
        [im_time, ~, ~, ~, im_chan, head] = ScanRead(name);
        %bin szie[s]
        resol = 5e-4;
        %number of time points
        timePoints = (head.TTResult_StopAfter/1000)/resol;
        %Discretize arrival times
        intensity = discretize(im_time,timePoints);
        %generate time trace
        timeTrace = accumarray(intensity,1);
        
        %PSD analysis
        [freq,PSD,degrees]=getPSD(timeTrace.',resol);
        
        plotPSD(freq,PSD,degrees,'b',ax)
        tmp = getFilename(name);
        tmp = append(tmp,'_',string(resol),'s');
        file_name = append(tmp,'_PSD.png');
        exportgraphics(ax, file_name,'Resolution',600)

    catch error
        disp(error.message)
    end
end


