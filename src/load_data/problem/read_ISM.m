function [im_chan,im_tcspc,im_posy,im_posx,im_time,head] = read_ISM(name)
%READ_ISM read 5d FLISM Data and header form PTU file 
%   Detailed explanation goes here
fname = name;

head  = PTU_Read_Head(fname);


%number of pixels in recorded image
s_pixl_x = head.ImgHdr_PixX;
s_pixl_y = head.ImgHdr_PixY;

time_R = mean(head.MeasDesc_Resolution);
head.max_bin = round(1./head.TTResult_SyncRate /time_R);

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
paras = [s_pixl_x s_pixl_y LStart LStop Frame head.ImgHdr_BiDirect head.length head.TTResultFormat_TTTRRecType(1)];

[im_data, im_time, im_param] = PTU_ReadPos(fname, paras);

%  im_data contains  detector-no, tcspc-bin, x-, and y-position
%  in 16 bit blocks
%  max number per block
block_max = 2^16-1;

im_data  = uint64(im_data);
im_chan  = double(bitand(bitshift(im_data,-48),block_max));
im_tcspc = double(bitand(bitshift(im_data,-32),block_max));
im_posy  = double(bitand(bitshift(im_data,-16),block_max)); % this is the line no.
im_posx  = double(bitand(im_data,block_max));
im_posx  = im_posx/block_max;

% im_param contains some statistics of the scan
head.ImgHdr_FrameNum  = im_param(1);  % how many frames were scanned
head.ImgHdr_LineNum   = im_param(2);  % how many lines were completed
head.ImgHdr_PixNum    = im_param(3);  % how many pixels wer scanned

head.ImgHdr_FrameTime = im_param(4)/head.TTResult_SyncRate;  % duration of one frame
head.ImgHdr_LineTime  = im_param(5)/head.TTResult_SyncRate;  % duration of one line
head.ImgHdr_PixelTime = im_param(6)/head.TTResult_SyncRate;  % accumulated time per pixel
head.ImgHdr_DwellTime = im_param(7)/head.TTResult_SyncRate;  % pixel dewll time during one scan
end

