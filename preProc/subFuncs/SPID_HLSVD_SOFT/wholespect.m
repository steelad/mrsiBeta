%***************************************************************************
%                                WHOLESPECT.M
%***************************************************************************
% PURPOSE:  Feature extraction by taking the whole spectrum in the region
% of interest
%***************************************************************************
% CALL:     [scores,misc] = wholespect(signal,step,frequency,ndp,begin,boundL,boundH)
%***************************************************************************
% INPUT:    signal      -- signal                        row vector
%           step        -- time step between points (ms) scalar
%           frequency   -- spectrometer frequency (kHz)  scalar
%           ndp         -- number of data points         scalar
%           begin       -- begin time (ms)               scalar
%           boundL      -- low bound                     scalar
%           boundH      -- high bound                    scalar
% OUTPUT:   scores      -- amplitude values of the points matrix
%           misc        -- see bottom of the file        structure
%**************************************************************************
function [scores,misc] = wholespect(signal,step,frequency,ndp,begin,boundL,boundH)
ppmaxis = ppmaxisfind(step,ndp,frequency);
ind = find((ppmaxis>=boundL)&(ppmaxis<=boundH)); 
[fftsig] = abs(fftshift(fft(signal)));
[scores] = fftsig(:,ind);
misc.nomisc = 1;