function [amp,tau] = double_exp(x,y)
%DOUBLE_EXP Fits double exponetial decay to x vs y data

max_amp = max(y);
init = [max_amp/2, max_amp/2, 2, 2, y(end)];
low = [0 0 0.1 0.1 0];
up = [max_amp max_amp 20 20 max_amp];
fitfun = fittype( @(amp1, amp2, tau1, tau2, offset,x) amp1*exp(-x/tau1) + amp2*exp(-x/tau2) + offset );
f = fit(x,y,fitfun,'StartPoint',init, 'Lower', low, 'Upper', up, 'Display', 'final');

amp = [f.amp1 f.amp2 f.offset];
tau = [f.tau1 f.tau2 inf];
end