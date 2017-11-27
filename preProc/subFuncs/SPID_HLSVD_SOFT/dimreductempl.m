%***************************************************************************
%                               dimreductempl.M
%***************************************************************************
% PURPOSE:  Feature selection using dimensionality reduction techniques
%***************************************************************************
% CALL:     [newscores,misc] = dimreductempl(signal,step,frequency,ndp,begin,...
procres, classopt,methname,nbrvar,boundL,boundH)
%***************************************************************************
% INPUT:    signal      -- signal                        row vector
%           step        -- time step between points (ms) scalar
%           frequency   -- spectrometer frequency (kHz)  scalar
%           ndp         -- number of data points         scalar
%           begin       -- begin time (ms)               scalar
%           procres     -- structure of processing results
%           classopt    -- structure of classification options
%           methname    -- method name                   string
%           nbrvar      -- number of variables to keep   integer
%           boundL      -- low bound                     scalar
%           boundH      -- high bound                    scalar
% OUTPUT:   newscores   -- new scores (best variables)   matrix
%           misc        -- see bottom of file            structure
%**************************************************************************
function [newscores,misc] = dimreductempl(signal,step,frequency,ndp,begin,procres,classopt,methname,nbrvar,boundL,boundH)
if ~isfield(classopt,'classtype')
    warning('The data are going to be mapped without taking into account the class they belong to.')
    %error('The class types are missing. ');
    return;
else
    if isempty(classopt.classtype)
        warning('The data are going to be mapped without taking into account the class they belong to.')
        %   error('The class types are missing. ');
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
        newscores(:,:,i) = compute_mapping(signal, methname, nbrvar);
        misc = 0;
    else
        compute_mapping(signal, methname);
        misc = 0;
    end
end


