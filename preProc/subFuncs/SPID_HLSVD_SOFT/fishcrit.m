%***************************************************************************
%                                fishcrit.M
%***************************************************************************
% PURPOSE:  Feature selection using the Fisher criterion
%***************************************************************************
% CALL:     [newscores,misc] = fishcrit(signal,procres,classopt,nbrvar)
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
%           misc(i).ranked_vector -- ranked vector       vector
%**************************************************************************
function [newscores,misc] = fishcrit(signal,step,frequency,ndp,begin,procres,classopt,nbrvar,boundL,boundH)

if nargin < 4
    nbrvar = 0;
end

if ~isfield(classopt,'classtype')
    error('Fisher criterion (supervised method) requires that a class has to be allocated to each data');
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
    group1 = scores(ind1,:);
    group2 = scores(ind2,:);

    ranked_vector = [];
    col_number = size(group1,2);
    temp = [];

    mean1 = [];
    mean2 = [];
    var1 = [];
    var2 = [];

    for ii = 1:col_number
        mean1 = [mean1 mean(group1(:,ii))];
        mean2 = [mean2 mean(group2(:,ii))];
        var1 = [var1 var(group1(:,ii))];
        var2 = [var2 var(group2(:,ii))];
    end
    for ii = 1:col_number
        val = ((abs(mean1(ii)-mean2(ii)))^2) / (var1(ii)+var2(ii));
        temp = [temp val];
    end

    [temp_sorted ranked_vector] = sort(temp,'descend');
    if nbrvar
        newscores(:,:,i) = scores(:,ranked_vector(1:nbrvar));
        misc(i).ranked_vector = ranked_vector(1:nbrvar);
    else
        newscores(:,:,i) = scores(:,ranked_vector);
        misc(i).ranked_vector = ranked_vector;
    end
end
