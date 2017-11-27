%***************************************************************************
%                                kruswal.M
%***************************************************************************
% PURPOSE:  Feature selection using the Kruskal-Wallis test
%***************************************************************************
% CALL:     [newscores,misc] = kruswal(signal,procres,classopt,nbrvar)
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
function [newscores,misc] = kruswal(signal,step,frequency,ndp,begin,procres,classopt,nbrvar,boundL,boundH)
if nargin < 4
    nbrvar = 0;
end

if ~isfield(classopt,'classtype')
    error('Fisher criterion (supervised method) requires that a class has to be allocated to each data');
    return;
else
    if isempty(classopt.classtype)
        error('''classtype'' unexisting: Fisher criterion (supervised method) requires that a class has to be allocated to each data');
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
    % ranked _vector contains in descending order the columns chi-square value ranking
    ranked_vector = [];
    temp = [];
    n_col = size(group1,2);

    for ii = 1:n_col
        [p,table] = kruskalwallis([group1(:,ii); group2(:,ii)],[-ones(size(group1(:,ii),1),1); ones(size(group2(:,ii),1),1)],'off');
        temp = [temp table(2,5)];
    end
    [temp_sorted ranked_vector] = sort(cell2mat(temp),'descend');
    if nbrvar
        newscores(:,:,i) = scores(:,ranked_vector(1:nbrvar));
        misc(i).ranked_vector = ranked_vector(1:nbrvar);
    else
        newscores(:,:,i) = scores(:,ranked_vector);
        misc(i).ranked_vector = ranked_vector;
    end
end
