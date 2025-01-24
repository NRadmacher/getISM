function [im_time, im_tcspc, im_posx, im_posy, im_chan, head] = ScanRead(name)
    head  = PTU_Read_Head(name);
    
    if head.TNTnPhoton == 0
        ME = MException('PTUread:nphoton','PTU file %s \ncontains no photons', name);
        throw(ME)
    end

    im_time  = [];
    im_tcspc = [];
    im_posx  = [];
    im_posy  = [];
    im_chan  = [];
    if head.ImgHdr_Dimensions == 1

        [im_sync, im_tcspc, im_chan, ~, ~, ~] = PTU_Read(name, head.TNTnPhoton, head);
        im_time = im_sync./head.TNTsyncRate;% + im_tcspc*head.TNTtcspsBinSize;
        return
    end
    nx     = head.ImgHdr_PixX;
    ny     = head.ImgHdr_PixY;
    
    LStart = 1;
    LStop  = 2;
    Frame  = 0;
    
    if isfield(head,'ImgHdr_LineStart')
        LStart = 2^(head.ImgHdr_LineStart-1);
    end
    if isfield(head,'ImgHdr_LineStop')
        LStop = 2^(head.ImgHdr_LineStop-1);
    end
    if isfield(head,'ImgHdr_Frame')
        Frame = 2^(head.ImgHdr_Frame-1);
    end
    
    [im_data, im_time, im_param] = PTU_ReadPos(name, [nx ny LStart LStop Frame head.ImgHdr_BiDirect head.length head.TTResultFormat_TTTRRecType(1)]);
    
    %  im_data contains  detector-no, tcspc-bin, x-, and y-position
    %  in 16 bit blocks
    
    im_data  = uint64(im_data);
    im_chan  = double(bitand(bitshift(im_data,-48),65535));
    im_tcspc = double(bitand(bitshift(im_data,-32),65535));
    im_posy  = double(bitand(bitshift(im_data,-16),65535)); % this is the line no.
    im_posx  = double(bitand(im_data,65535))/65536;         % this is the fraction of the current line.
    
    ind = (im_posy>0)&(im_posy<=ny);
    im_chan  = im_chan(ind);
    im_tcspc = im_tcspc(ind);
    im_posx  = im_posx(ind);
    im_posy  = im_posy(ind);            
    clear im_data;
        
    % im_param contains some statistics of the scan 
    head.ImgHdr_FrameNum  = im_param(1);  % how many frames were scanned
    head.ImgHdr_LineNum   = im_param(2);  % how many lines were completed
    head.ImgHdr_PixNum    = im_param(3);  % how many pixels wer scanned
    
    head.ImgHdr_FrameTime = im_param(4)/head.TTResult_SyncRate;  % duration of one frame
    head.ImgHdr_LineTime  = im_param(5)/head.TTResult_SyncRate;  % duration of one line
    head.ImgHdr_PixelTime = im_param(6)/head.TTResult_SyncRate;  % accumulated time per pixel
    head.ImgHdr_DwellTime = im_param(7)/head.TTResult_SyncRate;  % pixel dewll time during one scan
    
    time_R = mean(head.MeasDesc_Resolution);
    head.max_bin = round(1./head.TTResult_SyncRate /time_R);

end