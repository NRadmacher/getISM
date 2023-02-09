function plot_pixeldecay(tcspc,chan,n_pixl)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

max_bin = max(tcspc);
figure
hold on
for i = 1:n_pixl
    [c_count, ~] = histcounts(tcspc(chan == i-1), 1:max_bin+1);
    chr = int2str(i-1);
    plot(c_count, 'DisplayName',chr)
end
legend
end

