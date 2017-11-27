%***************************************************************************
%                                scalenorm.M
%***************************************************************************
% PURPOSE:  Normalize variables
%***************************************************************************
% CALL:     [scores,misc] = scalenorm(signal,step,frequency,ndp,begin,ir,freqmatrix,compNames)
%***************************************************************************
% INPUT:    signal      -- signal                        row vector
%           step        -- time step between points (ms) scalar
%           frequency   -- spectrometer frequency (kHz)  scalar
%           ndp         -- number of data points         scalar
%           begin       -- begin time (ms)               scalar
%           procres     -- structure of processing results
%           classopt    -- structure of classification options
%           scale       -- scaling option                integer
%           sigmoid     -- =1 : sigmoid shape            binary
% OUTPUT:   scores      -- spectrum point values         matrix
%           misc        -- see bottom of file            structure
%**************************************************************************
function [scores,misc] = scalenorm(signal, step, frequency,ndp,begin,procres,classopt,scale,sigmoid)

%  Hyperparameters, and their defaults
%  scale             -- 0: do nothing,
%                       1: columns have mean 0, std 1
%                       2: rows have mean, std 1
%  		        3: try to do do both 1 & 2
% !!!not available yet!!!  4: scale by correlation coefficients
%  sigmoid=0       -- yes or no, scale features nonlinearly through a
%                       sigmoid, i.e f(x)=tanh(sig*x)
%
if isempty(procres.scores)
    error('Features have to be extracted first and saved to procres.scores')
    return
end

x = procres.scores;
[l,n] = size(x);
%y = classopt.classtype;

if scale==4
    corr = zeros(size(x,2),1);
    corr = (mean(x(y==1,:))-mean(x(y==-1,:))).^2;
    st   = std(x(y==1,:)).^2+std(x(y==-1,:)).^2;
    a.corr = corr ./ st;
end

if scale==1 | scale==3  | scale==4
    a.mean_vec  = mean(x);
    a.scale_vec = std(x);
    if scale==3 %% try to do both scalings
        x = x - ones(l,1)*a.mean_vec; x = x * diag(1./a.scale_vec);
    end
end
if scale==2 | scale==3
    a.mean_vec2  = mean(x')';
    a.scale_vec2 = std(x')';
end


if scale==1 | scale==3 | scale==4
    x = x - ones(l,1)*a.mean_vec;
    x = x * diag(1./a.scale_vec);
end
if scale==2 | scale==3
    x = x - a.mean_vec2*ones(1,n);
    x = (x' * diag(1./a.scale_vec2))';
end

if sigmoid
    x=tanh(x*sigmoid);
end

if scale==4
    x = x * diag(a.corr);
end

scores = x;
misc = procres.misc;
