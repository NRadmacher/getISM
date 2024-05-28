function [head, im_tcspc, im_chan, im_line, im_col, im_time] = PTU_ScanRead(name, plt)

% if (nargin>1)&&~isempty(plt)
%     plt = 1;
% else
%     plt = 0;
% end

nphot = 5e7;
head  = [];

im_tcspc = [];
im_chan  = [];
im_line  = [];
im_col   = [];
im_time  = [];
im_frame = [];

if strcmp(name(end-2:end),'ptu')
    
    head = PTU_Read_Head(name);
    
    if ~isempty(head)
        
        nx = head.ImgHdr_PixX;
        ny = head.ImgHdr_PixY;
                
        if (head.ImgHdr_Ident == 1)||(head.ImgHdr_Ident == 6)
            
            anzch      = 32;
            Resolution = max(1e9*head.MeasDesc_Resolution);
            chDiv      = 1e-9*Resolution/head.MeasDesc_Resolution;
            Ngate      = ceil(1e9*head.MeasDesc_GlobalResolution./Resolution)+1;
            
            LineStart = 4;
            LineStop  = 2;
            
            if isfield(head,'ImgHdr_LineStart')
                LineStart = 2^(head.ImgHdr_LineStart-1);
            end
            if isfield(head,'ImgHdr_LineStop')
                LineStop = 2^(head.ImgHdr_LineStop-1);
            end
            
            y        = [];
            tmpx     = [];
            chan     = [];
            markers  = [];
            dt       = zeros(ny,1);
            
            Turns1   = [];
            Turns2   = [];
            
            cnt      = 0;
            tend     = 0;
            line     = 1;
            
            h = waitbar(0,'Please wait ...');
            
            fprintf('\n\n');
            
            if head.ImgHdr_BiDirect == 0
                
                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 nphot], head);
                
                while (num>0)
                    
                    cnt = cnt + num;
                    if ~isempty(y)
                        tmpy = tmpy+tend;
                    end
                    
                    ind = (tmpmarkers>0)|((tmpchan<anzch)&(tmptcspc<Ngate*chDiv));
                    
                    y       = [y; tmpy(ind)];                         %#ok<AGROW>
                    tmpx    = [tmpx; floor(tmptcspc(ind)./chDiv)+1;]; %#ok<AGROW>
                    chan    = [chan; tmpchan(ind)+1];                 %#ok<AGROW>
                    markers = [markers; tmpmarkers(ind)];             %#ok<AGROW>
                    
                    if LineStart==LineStop
                        tmpturns = y(markers==LineStart);
                        if numel(Turns1)>numel(Turns2)           % first turn is a LineStop
                            Turns1 = [Turns1; tmpturns(2:2:end)];      %#ok<AGROW>
                            Turns2 = [Turns2; tmpturns(1:2:end)];      %#ok<AGROW>
                        else
                            Turns1 = [Turns1; tmpturns(1:2:end)];      %#ok<AGROW>
                            Turns2 = [Turns2; tmpturns(2:2:end)];      %#ok<AGROW>
                        end
                    else
                        Turns1 = [Turns1; y(markers==LineStart)]; %#ok<AGROW>
                        Turns2 = [Turns2; y(markers==LineStop)];  %#ok<AGROW>
                    end
                    
                    ind          = (markers~=0);
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];
                    markers(ind) = [];
                    
                    tend  = y(end)+loc;
                    
                    if numel(Turns2)>1
                        for j=1:numel(Turns2)-1
                            
                            t1 = Turns1(1);
                            t2 = Turns2(1);
                            
                            ind          = (y<t1);
                            y(ind)       = [];
                            tmpx(ind)    = [];
                            chan(ind)    = [];
                            markers(ind) = [];
                            
                            ind = (y>=t1)&(y<=t2);
                            
                            im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              %#ok<AGROW>
                            im_chan   = [im_chan; uint8(chan(ind))];                %#ok<AGROW>
                            im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                            im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW>
                            im_time   = [im_time;  y(ind)];                          %#ok<AGROW>
                            
                            dt(line)  = t2-t1;
                            line = line +1;
                            waitbar(line/ny);
                            drawnow
                            
                            Turns1(1) = [];
                            Turns2(1) = [];
                        end
                    end
                    
                    [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 nphot], head);
                    
                end
                
                t1 = Turns1(end);
                t2 = Turns2(end);
                
                ind          = (y<t1);
                y(ind)       = [];
                tmpx(ind)    = [];
                chan(ind)    = [];
                
                ind = (y>=t1)&(y<=t2);
                
                im_tcspc  = [im_tcspc; uint16(tmpx(ind))];
                im_chan   = [im_chan; uint8(chan(ind))];
                im_line   = [im_line; uint16(line.*ones(sum(ind),1))];
                im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];
                im_time   = [im_time;  y(ind)];
                
                dt(line)  = t2-t1;
                
                line = line +1;
                waitbar(line/ny,h);
                drawnow
                
            else  % bidirectional scan
                
                rev = false;
                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 nphot], head);
                
                while (num>0)
                    
                    cnt = cnt + num;
                    if ~isempty(y)
                        tmpy = tmpy+tend;
                    end
                    
                    ind = ((tmpchan<anzch)&(tmptcspc<=Ngate*chDiv));
                    
                    y       = [y; tmpy(ind)];                         %#ok<AGROW>
                    tmpx    = [tmpx; floor(tmptcspc(ind)./chDiv)+1;]; %#ok<AGROW>
                    chan    = [chan; tmpchan(ind)+1];                   %#ok<AGROW>
                    markers = [markers; tmpmarkers(ind)];             %#ok<AGROW>
                    
                    if LineStart==LineStop
                        tmpturns = y(markers==LineStart);
                        if numel(Turns1)>numel(Turns2)           % first turn is a LineStop
                            Turns1 = [Turns1; tmpturns(2:2:end)];      %#ok<AGROW>
                            Turns2 = [Turns2; tmpturns(1:2:end)];      %#ok<AGROW>
                        else
                            Turns1 = [Turns1; tmpturns(1:2:end)];      %#ok<AGROW>
                            Turns2 = [Turns2; tmpturns(2:2:end)];      %#ok<AGROW>
                        end
                    else
                        Turns1 = [Turns1; y(markers==LineStart)]; %#ok<AGROW>
                        Turns2 = [Turns2; y(markers==LineStop)];  %#ok<AGROW>
                    end
                    
                    ind          = (markers~=0);
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];
                    markers(ind) = [];
                    
                    tend = y(end)+loc;
                    
                    for j=1:numel(Turns2)
                        
                        t1 = Turns1(1);
                        t2 = Turns2(1);
                        
                        ind          = (y<t1);
                        y(ind)       = [];
                        tmpx(ind)    = [];
                        chan(ind)    = [];
                        markers(ind) = [];
                        
                        ind = (y<=t2);
                        
                        im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              %#ok<AGROW>
                        im_chan   = [im_chan; uint8(chan(ind))];                %#ok<AGROW>
                        im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                        if rev
                            im_col    = [im_col;  uint16(nx - floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW>
                        else
                            im_col    = [im_col;  uint16( 1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW>
                        end
                        im_time   = [im_time;  y(ind)];                          %#ok<AGROW>
                        
                        dt(line)  = t2-t1;
                        
                        line = line +1;
                        rev  = ~rev;
                        
                        waitbar(line/ny,h);
                        drawnow
                        
                        Turns1(1) = [];
                        Turns2(1) = [];
                    end
                                        
                    [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 nphot], head);
                    
                end
                
                while ~isempty(Turns2)
                    t1 = Turns1(1);
                    t2 = Turns2(1);
                    
                    ind = (y<t1);        
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];
                    markers(ind) = [];
                    
                    ind = (y<=t2);
                    
                    im_tcspc  = [im_tcspc; uint16(tmpx(ind))];
                    im_chan   = [im_chan; uint8(chan(ind))];
                    im_line   = [im_line; uint16(line.*ones(sum(ind),1))];                   
                    if rev
                        im_col    = [im_col;  uint16(nx - floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW>
                    else
                        im_col    = [im_col;  uint16( 1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW>
                    end
                    im_time   = [im_time;  y(ind)];                          %#ok<AGROW>
                    
                    dt(line)  = t2-t1;
                    
                    line = line +1;
                    rev  = ~rev;
                    
                    Turns1(1) = [];
                    Turns2(1) = [];
                    
                    waitbar(line/ny);
                    drawnow
                end
            end

            im_frame  = ones(size(im_time));
            if ~isfield(head,'head.ImgHdr_PixelTime')
                head.ImgHdr_PixelTime = 1e9.*mean(dt(1:line-1))/nx/head.TTResult_SyncRate;
            end
            
            close(h);
            
        elseif (head.ImgHdr_Ident == 3)||(head.ImgHdr_Ident == 9)
            
            y        = [];
            tmpx     = [];
            chan     = [];
            marker   = [];
            
            dt       = zeros(ny,1);
            
            cnt      = 0;
            tend     = 0;
            line     = 1;
            n_frames = 0;
            f_times  = [];
            
            if ~isfield(head,'ImgHdr_X0')
                head.ImgHdr_X0 = 0;
            end
            if ~isfield(head,'ImgHdr_Y0')
                head.ImgHdr_Y0       = 0;
            end
            if ~isfield(head,'ImgHdr_PixResol')
                head.ImgHdr_PixResol = 1;
            end
            
            LineStart = 2^(head.ImgHdr_LineStart-1);
            LineStop  = 2^(head.ImgHdr_LineStop-1);
            Frame     = 2^(head.ImgHdr_Frame-1);
            
            h = waitbar(0, 'Please Wait...');
            
            if Frame < 1
                Frame = -1;
            end
            
            in_frame = true; %false;
            
            if Frame < 1
                in_frame = true;
                n_frames = n_frames + 1;
            end
            
            [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 nphot], head);
            
            while (num>0)
                
                t_tcspc  = [];
                t_chan   = [];
                t_line   = [];
                t_col    = [];
                t_y      = [];
                t_frame  = [];
                
                cnt = cnt + num;
                
                waitbar(cnt/head.TTResult_NumberOfRecords,h);
                
                tmpy = tmpy+tend;
                
                y       = [y; tmpy];                       %#ok<AGROW>
                tmpx    = [tmpx; tmptcspc];                %#ok<AGROW>
                chan    = [chan; tmpchan+1];               %#ok<AGROW>
                marker  = [marker; tmpmarkers];            %#ok<AGROW>
                tend    = y(end)+loc;
                
                F  = y(bitand(marker,Frame)>0);
                
                while ~isempty(F)
                    
                    if ~in_frame
                        ind = (y<=F(1));
                        y(ind)       = [];
                        tmpx(ind)    = [];
                        chan(ind)    = [];
                        marker(ind)  = [];
                        line         = 1;
                        in_frame     = true;
                        n_frames     = n_frames + 1;
                        f_times      = [f_times; F(1)];    %#ok<AGROW>
                        F(1)         = [];
                    end
                    
                    if ~isempty(F)
                        ind = y<F(1);
                        
                        f_y  = y(ind);
                        f_x  = tmpx(ind);
                        f_ch = chan(ind);
                        f_m  = marker(ind);
                        
                        y(ind)      = [];
                        tmpx(ind)   = [];
                        chan(ind)   = [];
                        marker(ind) = [];
                                                                        
                        L1 = f_y(bitand(f_m,LineStart)>0);
                        L2 = f_y(bitand(f_m,LineStop)>0);
                        
                        ll = line + numel(L2)-1; % this will be the last complete line in the data stack
                        
                        if ll > ny
                            L1 = L1(1:ny-line+1);
                            L2 = L2(1:ny-line+1);
                        end
                        
                        if numel(L1)>1
                            for j=1:numel(L2)
                                
                                ind = (f_y>L1(j))&(f_y<L2(j));
                                
                                t_tcspc  = [t_tcspc; uint16(f_x(ind))];              %#ok<AGROW>
                                t_chan   = [t_chan; uint8(f_ch(ind))];                %#ok<AGROW>
                                t_line   = [t_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                                t_col    = [t_col;  uint16(1 + floor(nx.*(f_y(ind)-L1(j))./(L2(j)-L1(j))))];  %#ok<AGROW>
                                t_y      = [t_y; f_y(ind)];                             %#ok<AGROW>
                                t_frame  = [t_frame;  uint16(n_frames.*ones(sum(ind),1))]; %#ok<AGROW>
                                dt(line) = dt(line) + (L2(j)-L1(j));
                                line = line +1;
                            end
                        end
                        
                        if line>ny
                            in_frame = false;
                        end
                    end                                               
                end
                
                im_tcspc  = [im_tcspc; t_tcspc];  %#ok<AGROW>
                im_chan   = [im_chan;  t_chan];   %#ok<AGROW>
                im_line   = [im_line;  t_line];   %#ok<AGROW>
                im_col    = [im_col;   t_col];    %#ok<AGROW>
                im_time   = [im_time;  t_y];      %#ok<AGROW>
                im_frame  = [im_frame; t_frame];      %#ok<AGROW>
                
                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = PTU_Read(name, [cnt+1 nphot], head);
                
            end
            
            F  = y(bitand(marker,Frame)>0);
            
            t_tcspc  = [];
            t_chan   = [];
            t_line   = [];
            t_col    = [];
            t_y      = [];
            t_frame  = [];
            
            if ~in_frame
                if isempty(F)
                    y       = [];
                    tmpx    = [];
                    chan    = [];
                    marker  = [];
                    line    = 1;
                else
                    ind = (y<=F(1));
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];
                    marker(ind)  = [];
                    line         = 1;
                    n_frames     = n_frames + 1;
                    f_times      = [f_times; F(1)];
                end
            end
            
            f_y = y;
            f_x = tmpx;
            f_ch = chan;
            f_m  = marker;
            
            clear y tmpx chan;
            
            L1 = f_y(bitand(f_m,LineStart)>0);
            L2 = f_y(bitand(f_m,LineStop)>0);
            
            ll = line + numel(L2)-1; % this will be the last complete line in the data stack
            if ll > ny
                L1 = L1(1:ny-line+1);
                L2 = L2(1:ny-line+1);
            end
            
            if numel(L1)>1
                for j=1:numel(L2)
                    ind = (f_y>L1(j))&(f_y<L2(j));
                    t_tcspc  = [t_tcspc; uint16(f_x(ind))];              %#ok<AGROW>
                    t_chan   = [t_chan; uint8(f_ch(ind))];                %#ok<AGROW>
                    t_line   = [t_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                    t_col    = [t_col;  uint16(1 + floor(nx.*(f_y(ind)-L1(j))./(L2(j)-L1(j))))];  %#ok<AGROW>
                    t_y      = [t_y;    f_y(ind)];  %#ok<AGROW>
                    t_frame  = [t_frame;  uint16(n_frames.*ones(sum(ind),1))]; %#ok<AGROW>
                    dt(line) = dt(line) + (L2(j)-L1(j));
                    line = line +1;
                end
            end
            
            im_tcspc  = [im_tcspc; t_tcspc];
            im_chan   = [im_chan;  t_chan];
            im_line   = [im_line;  t_line];
            im_col    = [im_col;   t_col];
            im_time   = [im_time;  t_y];
            im_frame  = [im_frame; t_frame];      
            
            if head.ImgHdr_BiDirect == 1
                ind = (mod(im_line,2)==0);
                im_col(ind) = 1 + nx-im_col(ind);
            end 
            
            head.ImgHdr_FrameTime = 1e9.*mean(diff(f_times))/head.TTResult_SyncRate;
            head.ImgHdr_PixelTime = 1e9.*mean(dt)/nx/head.TTResult_SyncRate;
            head.ImgHdr_DwellTime = head.ImgHdr_PixelTime./n_frames;
            
            close(h);
        end
%        save([name(1:end-4) '_raw.mat'],'im_tcspc','im_chan','im_line','im_col','im_time','im_frame','-v7.3');
    end
end

im_time = im_time./head.TTResult_SyncRate + double(im_tcspc).*head.MeasDesc_Resolution; % in seconds
%~isempty(head)&&(plt==1)
if plt==1

    xl = 1:255/7:256;
    yl = [[0 0 0]; [0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]; [1 0 1]];
    lambda   = 1:256;
    spectrum = interp1(xl, yl, lambda);

    tag = zeros(ny,nx);
    for y = 1:ny
        ind = (im_line == y);
        tmp1 = im_col(ind);
        for x = 1:nx
            tag(y,x) = sum(tmp1 == x);
        end
    end
    
    figure
    caxis([min(tag(:)) max(tag(:))])
    colormap(spectrum)
    
    x = head.ImgHdr_X0+(1:nx)*head.ImgHdr_PixResol;
    y = head.ImgHdr_Y0+(1:ny)*head.ImgHdr_PixResol;
    imagesc(x,y,tag);
    set(gca,'DataAspectRatio', [1,1,1], ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'XDir','normal', ...
        'YDir','reverse');
    xlabel(sprintf('x / %sm',char(181)));
    ylabel(sprintf('y / %sm',char(181)));
end