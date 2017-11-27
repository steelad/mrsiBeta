%***************************************************************************
%                                peakbuckettempl.M
%***************************************************************************
% PURPOSE:  Feature extraction with peak bucketing 
%***************************************************************************
% CALL:     [scores,misc] = peakbuckettempl(signal,step,frequency,ndp,begin,...
% ir,freqmatrix,compNames)
%***************************************************************************
% INPUT:    signal      -- signal                        row vector
%           step        -- time step between points (ms) scalar
%           frequency   -- spectrometer frequency (kHz)  scalar
%           ndp         -- number of data points         scalar
%           begin       -- begin time (ms)               scalar
%           ir          -- interval types                scalar (integer)
%           freqmatrix  -- interval matrix               matrix
%           compNames   -- Names of the intervals        cell
% OUTPUT:   scores      -- peak integrated values        matrix
%           misc        -- see bottom of the file        structure
%**************************************************************************
function [scores,misc] = peakbuckettempl(signal,step,frequency,ndp,begin,ir,freqmatrix,compNames)
if nargin < 7
    freqmatrix = [];
else
    if ~isempty(freqmatrix)
        for i = 1:length(freqmatrix)/2
            freqmatrixt(i,1:2) = freqmatrix((i-1)*2+1:i*2);
        end
        freqmatrix = freqmatrixt;
    end
end
if  nargin < 8
    compNames = {};
end
if (ir == 0) |(ir == 5) % user's intervals
    if isempty(freqmatrix)
        error('The intervals where the signals have to be integrated have to be defined')
    end
    if isempty(compNames)
        error('The list of components has to be defined')
    end
else
    [compNames,freqmatrix] = readintervals(ir);
end
ppmaxis = ppmaxisfind(step,ndp,frequency);
scores = [];
fftsig = abs(fftshift(fft(signal)));
spacing = ppmaxis(2)-ppmaxis(1);
for m = 1:size(freqmatrix,1)
    ind = find(ppmaxis-(spacing/2)<=(freqmatrix(m,2)) & ppmaxis+(spacing/2)>=(freqmatrix(m,1)));
    scores = [scores fftsig(ind)];
end
misc.compNames = compNames;
misc.peaks = freqmatrix;
