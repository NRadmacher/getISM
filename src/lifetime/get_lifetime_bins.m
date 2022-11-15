function [tail_t, tail_bin_l, tail_start_time, tcspc_t] = get_lifetime_bins(tail_start, tail_end, time_R, bin_factor)
%GET_LIFETIME_BINS retuns tcspc bin values for lifetime tail fitting for

% TCSPC bin length in ns
tcspc_bin_l = time_R * 1e9;

% TCSPC tail start in ns
tail_start_time = tail_start * tcspc_bin_l;

% tail edeges after binning in TCSPC indecis
tcspc_t     = tail_start:bin_factor:tail_end;

% tail length and bins in ns for tail fit
tail_bin_l  = tcspc_bin_l*bin_factor;
tail_t      = tcspc_bin_l*tcspc_t(1:end-1);
end