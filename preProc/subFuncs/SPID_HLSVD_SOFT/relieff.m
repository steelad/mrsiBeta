%***************************************************************************
%                                relieff.M
%***************************************************************************
% PURPOSE:  Feature selection using  Relief-f
%***************************************************************************
% CALL:     [newscores,misc] = relieff(signal,procres,classopt,nbrvar,iterations,neighbors)
%***************************************************************************
% INPUT:    signal      -- signal                        row vector
%           step        -- time step between points (ms) scalar
%           frequency   -- spectrometer frequency (kHz)  scalar
%           ndp         -- number of data points         scalar
%           begin       -- begin time (ms)               scalar
%           procres     -- structure of processing results
%           classopt    -- structure of classification options
%           nbrvar      -- number of variables to keep   integer
%           iterations  -- number of iterations          scalar
%           neighbors   -- number of neighbors           scalar   
%           boundL      -- low bound                     scalar
%           boundH      -- high bound                    scalar
% OUTPUT:   newscores   -- new scores (best variables)   matrix
%           misc(i).ranked_vector -- ranked vector       vector
%**************************************************************************
function [newscores,misc] = relieff(signal,step,frequency,ndp,begin,procres,classopt,nbrvar,iterations,neighbors,boundL,boundH)
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
    temp = [];
    total_data = [group1;group2;];
    row_number = size(total_data,1);
    row_number1 = size(group1,1);
    row_number2 = size(group2,1);
    col_number = size(group1,2);
    weights = zeros(1,col_number);
    max_values = max(total_data);
    min_values = min(total_data);

    if (min([row_number1 row_number2]) <= neighbors)
        neighbors = min([row_number1 row_number2]) - 1;
    end

    for ii = 1:iterations
        numbers = randperm(row_number);
        current = numbers(1);
        [hits,missess] = find_neighbors(group1,group2,current,neighbors);
        if (current > size(group1,1))
            option = 2;
        else
            option = 1;
        end;
        for a = 1:col_number

            if option == 1
                term_hits = diff_sum(total_data(current,:),group1,hits,a,max_values,min_values);
                term_missess = diff_sum(total_data(current,:),group2,missess,a,max_values,min_values);
            else
                term_hits = diff_sum(total_data(current,:),group2,hits,a,max_values,min_values);
                term_missess = diff_sum(total_data(current,:),group1,missess,a,max_values,min_values);
            end
            weights(a) = weights(a) - (term_hits/(iterations*neighbors)) + (term_missess/(iterations*neighbors));
        end
    end

    [weights_sorted ranked_vector] = sort(weights,'descend');
    if nbrvar
        newscores(:,:,i) = scores(:,ranked_vector(1:nbrvar));
        misc(i).ranked_vector = ranked_vector(1:nbrvar);
    else
        newscores(:,:,i) = scores(:,ranked_vector);
        misc(i).ranked_vector = ranked_vector;
    end
end
