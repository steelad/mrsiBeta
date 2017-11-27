%***************************************************************************
%                               ICAtemp.M
%***************************************************************************
% PURPOSE:  Feature selection using ICA
%***************************************************************************
% CALL:     [newscores,misc] = ICAtempl(signal,procres,classopt,nbrvar)
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
function [newscores,misc] = ICAtempl(signal,step,frequency,ndp,begin,procres,classopt,nbrvar,boundL,boundH)
if ~isfield(classopt,'classtype')
    error('ICA requires that a class is allocated to each data. The selection of the variables depends on the pair of classes.');
    return;
else
    if isempty(classopt.classtype)
        error('ICA requires that a class has to be allocated to each data');
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
    if nbrvar
        [zscores, A, W] = fastica(scores,'numOfIC', nbrvar,'displayMode','off','verbose','off');
    else
        [zscores, A, W] = fastica(scores,'displayMode','off','verbose','off');
    end
    newscores(:,:,i) = scores*zscores';
    misc.estsepmat = W;
    misc.mixmat = A;
end
