function [demptype,gausstype,freqtype,etype] = ...
    priorknowledge(nrstofs,prior_equal,equal_to)
% minimal priorknowledge

% initialize
demptype  = zeros(nrstofs,3);
gausstype = zeros(nrstofs,3);
freqtype  = zeros(nrstofs,3);
etype     = zeros(nrstofs,3);

% set:
% some equal dempings
demptype(prior_equal,1) = ones(length(prior_equal),1);
demptype(prior_equal,2) = equal_to';
% some equal frequencies
freqtype(prior_equal,1) = ones(length(prior_equal),1);
freqtype(prior_equal,2) = equal_to';
