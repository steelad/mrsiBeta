function W = ...
    init_weightingMRSI(demp,nra,nrd,nrg,nrf,nrph,nre,...
                 loline,galine,voline,weightfordamp,...
                 weightforgauss,weightforfreq,sigma);
% function that initializes the weighting matrix W needed for
% specifying weights of different types of parameters used
% in MRSI penalty term, through the parameters weightfordamp,
% weightforgauss,weightforfreq
% sigma is the noise standard deviation, included in W for convenience
%

nrstofs  = length(demp);

% initialize zero weights for amplitude variables
Wa = zeros(nra,1);

% initialize weights for demping variables
Wd = [];
if (loline | voline)
    Wd = weightfordamp*ones(nrd,1);
end

% initialize weights for gauss line variables
Wg = [];
if (galine | voline)
    Wg = weightforgauss*ones(nrg,1);
end

% initialize weights for frequency variables
Wf = weightforfreq*ones(nrf,1);

% initialize zero weights for phase variable
Wph = zeros(nrph,1);

% initialize zero weights for eddy current variables
We = zeros(nre,1);

W = sigma*diag([Wa;Wd;Wg;Wf;Wph;We]);
