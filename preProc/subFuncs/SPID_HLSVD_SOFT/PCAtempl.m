%***************************************************************************
%                               PCAtemp.M
%***************************************************************************
% PURPOSE:  Feature selection using PCA
%***************************************************************************
% CALL:     [newscores,misc] = PCAtempl(signal,procres,classopt,nbrvar)
%***************************************************************************
% INPUT:    signal      -- signal                        row vector
%           step        -- time step between points (ms) scalar
%           frequency   -- spectrometer frequency (kHz)  scalar
%           ndp         -- number of data points         scalar
%           begin       -- begin time (ms)               scalar
%           procres     -- structure of processing results
%           classopt    -- structure of classification options
%           nbrvar      -- number of variables to keep   integer
%           boundL      -- low bound                     scalar
%           boundH      -- high bound                    scalar
% OUTPUT:   newscores   -- new scores (best variables)   matrix
%           misc        -- see bottom of file            structure
%**************************************************************************
function [newscores,misc] = PCAtempl(signal,step,frequency,ndp,begin,procres,classopt,nbrvar,boundL,boundH)
if ~isfield(classopt,'classtype')
    error('PCA requires that a class is allocated to each data. The selection of the variables depends on the pair of classes.');
    return;
else
    if isempty(classopt.classtype)
        error('Fisher criterion (supervised method) requires that a class has to be allocated to each data');
        return;
    else
        classtype = classopt.classtype;
    end
end

[k,comb] = binaryencoding(classtype);
for i = 1:comb
    ind1 = find(k(:,i) == 1);
    ind2 = find(k(:,i) == 2);
    if isfield(procres,'scores')
        if length(size(procres.scores))==3
            scores = procres.scores(:,:,i);
        else
            scores = procres.scores;
        end
    else
        [xpts,xpts2] = absi(ndp,step,frequency);
        id = find(xpts2>boundL&xpts2<boundH);
        scores = signal(:,id);
    end
    [pc, zscores, pcvars] = princomp(scores,'econ');
    if nbrvar
        newscores(:,:,i) = zscores(:,1:nbrvar);
        misc(i).pc = pc(:,1:nbrvar);
    else
        newscores(:,:,i) = zscores;
        misc(i).pc = pc;
    end
end
