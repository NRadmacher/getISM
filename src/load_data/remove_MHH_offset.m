function [new_tcspc] = remove_MHH_offset(im_tcspc,im_chan,max_bin, hm_bin)
%REMOVE_MHH_OFFSET Removes time offset between two eventtimers
%   estimated time of half max, rising edge of tcspc curve
%   max_bin hightest tcspc bin
%   hm_bin  new pos of rising half maximum
m_time = zeros(23,1);

for i = 1:23
    [c_count, ~] = histcounts(im_tcspc(im_chan == i-1), 1:max_bin);
    % Find the half max value.
    halfMax = (min(c_count) + max(c_count)) / 2;
    m_time(i) = find(c_count >= halfMax, 1, 'first' );
end

shift_time = m_time - hm_bin;

for i = 1:23
    ind = im_chan == i-1;
    im_tcspc(ind) = mod(im_tcspc(ind) - shift_time(i), max_bin);
end
new_tcspc = im_tcspc;
end

