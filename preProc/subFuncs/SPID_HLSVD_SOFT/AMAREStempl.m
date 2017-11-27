%***************************************************************************
%                                AMAREStempl.M
%***************************************************************************
% PURPOSE:  Quantitation using AMARES
%***************************************************************************
% CALL:     [scores,misc] = AMAREStempl(signal,step,frequency,ndp,begin,MPFIR)
%***************************************************************************
% INPUT:    signal      -- signal                        row vector
%           step        -- time step between points (ms) scalar
%           frequency   -- spectrometer frequency (kHz)  scalar
%           ndp         -- number of data points         scalar
%           begin       -- begin time (ms)               scalar
%           MPFIR       -- MPFIR =1, nothing =0          binary
% OUTPUT:   newscores   -- new scores (best variables)   matrix
%           misc        -- nothing                       structure
%**************************************************************************
function [scores,misc] = AMAREStempl(signal,step,frequency,ndp,begin,MPFIR)
save sim signal step frequency ndp begin
spiddir = getappdata(0, 'spiddir');
prog_dir=fullfile(spiddir,'Software','processing');
eval(['!cd ',prog_dir])
if MPFIR % AMARES with MP-FIR (method by Sundin)
    %pg = fullfile(prog_dir,'amaresfreq');
    pg = fullfile(prog_dir,'filteramares');
else % AMARES without MP-FIR (method by Sundin)
    pg = fullfile(prog_dir,'amares');
end
eval(['!',pg,' sim.mat']);  %the signal file (.mat) has to be specified if it's not contained in the paramater file (.par)
%load lorentz
%eval(['!rm sim.mat'])
%eval(['!rm lorentz'])
misc = 0;
scores =0;

