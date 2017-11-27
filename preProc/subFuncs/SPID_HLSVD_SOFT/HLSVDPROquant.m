%***************************************************************************
%                                HLSVDPROquant.M
%***************************************************************************
% PURPOSE:  Filter out a specific region of the spectrum using the
% HLSVD-PRO by T. Laudadio et al.
%***************************************************************************
% CALL:     signal = HLSVDPROquant(signal,step,frequency,ndp,begin,boundL,boundH,M)
%***************************************************************************
% INPUT:    signal      -- signal                        row vector
%           step        -- time step between points (ms) scalar
%           frequency   -- spectrometer frequency (kHz)  scalar
%           ndp         -- number of data points         scalar
%           begin       -- begin time (ms)               scalar
%           boundL      -- low bound                     scalar
%           boundH      -- high bound                    scalar
%           M           -- model order                   scalar
% OUTPUT:   scores      -- filtered signal               row vector
%           misc        -- see bottom of the file        structure
%**************************************************************************
function [scores,misc] = HLSVDPROquant(signal,step,frequency,ndp,begin,boundL,boundH,M)

t = [0:step:(ndp-1)*step];
fsampl = 1/step;
xmin = (boundL-4.7)*frequency/10^6;%in kHz
xmax = (boundH-4.7)*frequency/10^6;%in kHz

freqrange1 = [-(fsampl*1000/2-1),xmin*1000]; %in Hz
freqrange2 = [xmax*1000,(fsampl*1000/2-1)]; %in Hz
freqrange = [freqrange1;freqrange2];
MM = round(ndp/2);
fsampl = 1/step*1000; %in Hz
[signal,freq,damp,ampl,phas,SV] = HLSVDPROrec(signal,M,fsampl,t/1000,MM,freqrange);
misc.ampl = ampl; %amplitude estimates 
misc.freq = freq; %frequency estimates
misc.phas = phas*180/pi; %phase estimates in degrees
misc.damp = damp; %damping estimates
misc.SV = SV; %singular values 
misc.filrecon=reconsd(t/1000,misc.freq,misc.damp,misc.ampl,misc.phas); %reconstructed signals
misc.filoriginal = signal; % these signals are the original ones   
for k = 1:M
      individual_components(k,:) = reconsd(t/1000,freq(k),damp(k),ampl(k),phas(k));
end
misc.individual_components = individual_components;
misc.baseline = [];
misc.residual =  misc.filrecon-misc.filoriginal;
misc.filter_information.order = M;
scores  = ampl;
