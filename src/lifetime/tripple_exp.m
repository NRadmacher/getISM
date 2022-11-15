function [amp,tau] = tripple_exp(x,y,plt)
%UNTITLED2 Summary of this function goes here

max_amp = max(y);
init = [max_amp/3, max_amp/3, max_amp/3, 0.5, 2, 4, 20];
low = [0 0 0  0 0 0];
up = [max_amp max_amp max_amp 20 20 20 max_amp];
fitfun = fittype( @(amp1, amp2, amp3, tau1, tau2, tau3, offset,x) amp1*exp(-x/tau1) + amp2*exp(-x/tau2) + amp3*exp(-x/tau3) + offset );
f = fit(x,y,fitfun,'StartPoint',init, 'Lower', low, 'Upper', up, 'Display', 'final');

if plt
    figure
    plot(f, x,y)
    set(gca, 'YScale', 'log')
    xlabel('time [ns]')
    ylabel('normalized count')
    
    tau_text = compose('tau: %0.1f %0.1f %0.1f', f.tau1, f.tau2, f.tau3);
    amp_text = compose('amp: %0.1f %0.1f %0.1f', f.amp1, f.amp2, f.amp3);
    
    text(max(x)/2,0.9*max(y),tau_text);
    % text(max(x)/2,0.7*max(y),amp_text);
end

amp = [f.amp1 f.amp2 f.amp3 f.offset];
tau = [f.tau1 f.tau2 f.tau3];

end