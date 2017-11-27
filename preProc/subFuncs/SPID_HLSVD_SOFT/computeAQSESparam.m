function [allow_damp,allow_freq,lambda,maxiter,phasedistort,filterlength,ripple,baseline_boolean,plotting] = computeAQSESparam(step,begin,ndp,frequency)
lambda = 2/step;
allow_damp = step/10;
allow_freq = step/10;
filterlength = 50;
maxiter = 25;
phasedistort = 0; % 0 => does not take into account the eddy current effect 
ripple = 0.01;
baseline_boolean = 0;
plotting = 0;


