function [alf,ampl,demp,gauss,freq,phas,e,lb,ub,nra,nrd,nrg,nrf,nrph,nre,...
          amplref,dempref,gaussref,freqref,phasref,eref] = ...
    init_one_run(ampl,demp,gauss,freq,phas,e,...
                 demptype,gausstype,freqtype,etype,...
                 loline,galine,voline,phasedistort,equalph,...
                 allow_change_damp,allow_change_freq);
% function that initializes the parameter vector "alf" and the
% lower and upper bounds for all variables
%
% In the case equalph = 0,
%    only the nonlinear variables (damping and frequency) are counted 
%    due to the VARPRO method used as NLLS algor.
% else, the amplitudes and the equal phase are also initialized.
%
%
% NOTE: allow_change_* defines the percentage in which 
%       the demp and gauss are allowed to change


nrstofs  = length(demp);

amplref  = [];
dempref  = [];
gaussref = [];
freqref  = [];
phasref  = [];
eref     = [];

index = 0;

lb = zeros(6*nrstofs,1);
ub = zeros(6*nrstofs,1);

% initialize amplitude variables
nra = 0;
if (equalph)
    for i = 1:nrstofs
      nra         = nra+1;
      index       = index+1;
      amplref(i)  = index;
      alf(index)  = ampl(i);
      temp(index) = 1;
      lb(index)   = 0; 
      ub(index)   = inf;
    end
end

% initialize demping variables
if (length(allow_change_damp)==1) 
    allow_change_damp = allow_change_damp(ones(2*nrstofs,1)); 
end
nrd = 0;
if (loline | voline)
    for i = 1:nrstofs
        select	= demptype(i,1);
        if (select == 0)               % free value for a demping
            nrd         = nrd+1;
            index       = index+1;
            dempref(i)  = index;
            alf(index)  = demp(i);
            temp(index) = 1;
            lb(index)   = max(0,-allow_change_damp(i)); 
            ub(index)   = allow_change_damp(nrstofs+i);
        elseif (select == -2)          % fixed value for a demping
            demp(i)     = demptype(i,3);
    end
end
for i = 1:nrstofs
    select = demptype(i,1);
    if (select == 1)                   % shift value for a demping
            ref         = dempref(demptype(i,2));
            dempref(i)  = ref;
            temp(ref)   = temp(ref)+1;
            alf(ref)    = alf(ref) + demp(i) - demptype(i,3);
        end
    end
end

% initialize gauss line variables
nrg = 0;
if (galine | voline)
    select  = gausstype(i,1);
    for i = 1:nrstofs
        select = gausstype(i,1);
        if (select == 0)                % free variable 
            nrg         = nrg+1;
            index       = index+1;
            gaussref(i) = index;
            alf(index)  = gauss(i);
            temp(index) = 1;
            lb(index)   = -inf; 
            ub(index)   = inf;
        elseif (select == -2)           % fixed
            gauss(i)    = gausstype(i,3);
    end
end
for i = 1:nrstofs
    select = gausstype(i,1);
    if (select == 1)                    % shift
            ref         = gaussref(gausstype(i,2));
            gaussref(i) = ref;
            temp(ref)   = temp(ref)+1;
            alf(ref)    = alf(ref) + gauss(i) - gausstype(i,3);
        end
    end
end



% initialize frequency variables
if (length(allow_change_freq)==1) 
    allow_change_freq = allow_change_freq(ones(2*nrstofs,1)); 
end
nrf = 0;
for i = 1:nrstofs
    select  = freqtype(i,1);
    if (select == 0)                    % free variable
      nrf         = nrf+1;
      index       = index+1;
      freqref(i)  = index;
      alf(index)  = freq(i);
      temp(index) = 1;
      lb(index) = -allow_change_freq(i);
      ub(index) = allow_change_freq(nrstofs+i);
    elseif (select == -2)               % fixed          
      freq(i)     = freqtype(i,3);
    end
end
for i = 1:nrstofs
    select = freqtype(i,1);
    if (select == 1)                    % shift
      ref         = freqref(freqtype(i,2));
      freqref(i)  = ref;
      temp(ref)   = temp(ref)+1;
      alf(ref)    = alf(ref) + freq(i) - freqtype(i,3);
    end
end

% initialize phase variable
nrph = 0;
if (equalph)
  nrph        = 1;
  index       = index+1;
  for i = 1:nrstofs
    phasref(i)  = index;
  end
  [maxampl,maxamplind] = max(ampl);
  alf(index)  = phas(maxamplind);
  temp(index) = 1;
  lb(index)   = -pi; 
  ub(index)   = pi;
end

% initialize eddy current variables
nre = 0;
if (phasedistort)
    for i = 1:nrstofs
        select = etype(i,1);
        if (select == 0)                 % free variable
            nre         = nre+1;
            index       = index+1;
            eref(i)     = index;
            alf(index)  = e(i);
            temp(index) = 1;
            lb(index)   = -inf;
            ub(index)   = +inf;
        elseif (select  == -2)           % fixed  
            e(i)        = etype(i,3);
        end
    end
    for i = 1:nrstofs
        select = etype(i,1);
        if (select == 1)                 % shift
            ref       = eref(etype(i,2));
            eref(i)   = ref;
            temp(ref) = temp(ref)+1;
            alf(ref)  = alf(ref) + e(i) - etype(i,3);
        end
    end
end

ub(index+1:end) = [];
lb(index+1:end) = [];

alf=alf./temp;