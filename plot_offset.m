%% 
clear
close all
clc

%% Loadind data

fname = 'C:\Users\NRadmacher\Documents\Uni\PHD\Messung_Daten\210909\offset_012.ptu';

[im_chan,im_tcspc,~,~,~, head] = read_ISM(fname);

m_time = zeros(23,1);
figure
hold on
for i = 1:23
    [c_count, time] = histcounts(im_tcspc(im_chan == i-1), 1:1275);
    m_time(i) = find(c_count >= max(c_count)/2, 1 );
%     FX = gradient(c_count);
    chr = int2str(i-1);
%     legend(chr)
    plot(c_count, 'DisplayName',chr)
%     plot(m_time(i), c_count(m_time(i)), 'o')
end
legend
m_time = m_time(m_time~=1);
mean_one = mean(m_time(1:15));
mean_two = mean(m_time(17:end));
shift_time = mean_one - mean_two;

for i = 17:23
    ind = im_chan == i-1;
    im_tcspc(ind) = im_tcspc(ind) + shift_time;
    
end

figure
hold on
for i = 1:23
    [c_count, ~] = histcounts(im_tcspc(im_chan == i-1), 1:1275);
    chr = int2str(i-1);
    plot(c_count, 'DisplayName',chr)
end
legend

max_bin = max(max(im_tcspc));

[count, ~] = histcounts(im_tcspc, 1:max_bin+1);
count = count(1:end-100);

figure
hold on
plot(count, '.')

figure
plot(m_time, 'x')
