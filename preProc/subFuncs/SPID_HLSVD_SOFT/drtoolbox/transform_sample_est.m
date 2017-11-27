function t_point = transform_sample_est(point, X, mappedX)
%TRANSFORM_SAMPLE_EST Performs out-of-sample extension using estimation technique
%
%   t_point = transform_sample_est(point, X, mappedX)
%
% Performs out-of-sample extension using estimation technique on datapoint
% points. You also need to specify the original dataset in X, and the
% reduced dataset in mappedX (the two datasets may also be PRTools datasets).
% The function returns the coordinates of the transformed point in t_point.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.1b.
% The toolbox can be obtained from http://www.cs.unimaas.nl/l.vandermaaten
% You are free to use, change, or redistribute this code in any way you
% want. However, it is appreciated if you maintain the name of the original
% author.
%
% (C) Laurens van der Maaten
% Maastricht University, 2007

    if strcmp(class(X), 'dataset')
        X = X.data;
    end
    if strcmp(class(mappedX), 'dataset')
        mappedX = mappedX.data;
    end

    % Find nearest neighbor for point
    n = size(X, 1);
    D = zeros(1, n);
    aa = sum(point .* point);
    bb = sum(X' .* X');
    ab = point * X';
    d = sqrt(repmat(aa', [1 size(bb, 2)]) + repmat(bb, [size(aa, 2) 1]) - 2 * ab);
    [d, ind] = min(d);

    % Compute transformation matrix
    L = pinv(X(ind,:) - mean(X(ind,:))) * (mappedX(ind,:) - mean(mappedX(ind,:)));
    
    % Compute coordinates of transformed point
    t_point = (mean(mappedX(ind,:)) + ((point - mean(X(ind,:))) * L))';
    