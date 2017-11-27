function [recon,baseline,filrecon,filbaseline,original,...
          filter_information,ampl,demp,phas,freq,gauss,e,...
          ampl_ratios,ampl_absolute,CR_error,individual_components, ...
          residual] = ...
    akses_gcv(invivosignal,dbase,t,step,abs_scl_factor,...
              demp0,freq0,gauss0,e0,loline,galine,voline, ...
              phasedistort,filtertype,flow,fhigh,ripple,filterlength,...
              prior_equal,equal_to,allow_change,baseline_boolean,lambda,...
              maxiter,plotting)
  
% Function that fits a short echo time proton spectrum in the time domain 
% using a metabolite database.
%
% Method: * nonlinear least squares with prior knowledge (linear
%           relations between some parameters) and linear bounds.
%         * automatic choice of the regularization parameter by
%           Generalized Cross Validation
%
% Input variables:
%   invivosignal  = signal to be quantified
%   dbase         = data base of metabolite signals, stored in columns
%   t             = vector of time instances 
%                   (default [0:length(invivosignal)-1])
%   step          = frequency sampling step (default 1)
%   abs_scl_facto = scaling factors for amplitudes (default ones)
%   demp0, freq0, gauss0, e0 = initial guesses for nonlinear parameters
%                              (default zeros(size(dbase,2),1) 
%   loline, galine, voline   = 0 or 1, to specify the line shapes - lorentz, 
%                              gauss or voigt (default 1, 0, 0)
%   phasedistort  = 0 or 1, to specify whether Eddy current correction is done
%   filtertype    = 'fir', for FIR filter, 'hlsvd', for HLSVDPRO,
%                    otherwise, no filter
%   flow,fhigh    = passband limits for filtering (in kHz)
%   ripple        = ripple suppresion when using a FIR filter (default 0.01)
%   filterlength  = minimal filterlength when using a FIR filter 
%   prior_equal, equal_to = vectors of indeces specifying simple priorknowledge
%                           i.e., dampings and frequences of metabolites
%                           in prior_equal are equal to equal_to
%                           (default [], [])
%   allow_change  = constant value that fixes the linear bounds on
%                   dampings and frequencies (default 0.1)
%   baseline_boolean = flag: 0 - no baseline is fit (default) 
%                            1 - baseline is fit
%   lambda        = regularization parameter for baseline coefficients
%                   (default 0.01)
%   maxiter       = maximum number of Levenberg-Marquardt iterations
%                   (default 10)
%   plotting      = flag: 0 - no plots during optimization process
%                         1 - plot frequency domain fit at every iteration
%
% Output variables:
%   recon         = reconstructed (fitted) signal
%                   length:  length(invivosignal);
%   baseline      = fitted baseline (time-domain)
%                   length:  length(invivosignal);
%   filrecon      = filtered reconstructed (fitted) signal
%                   length:  length(invivosignal) - length(fir_h) + 1;
%   filbaseline   = filtered fitted baseline (time-domain)
%                   length:  length(invivosignal) - length(fir_h) + 1;
%   original      = filtered original in vivo signal 
%                   length:  length(invivosignal) - length(fir_h) + 1;
%   ampl, demp, phas, freq, gauss, e = parameters of fitted signal
%                   length:  nrbasis = size(dbase,2);
%   ampl_ratios   = amplitude ratios to Cr amplitude
%   ampl_absolute = metabolite concentrations
%   CR_error      = Cramer-Rao error bounds 
%                   a structure containing the fields: a, d, g, f, ph, e,
%                   such that, e.g., CR_error.a is a vector of CR bounds
%                   for the amplitudes
%   individual_components = filtered individual components in data base
%   residual      = residual of the filtered fit


% Diana Sima, KUL-ESAT, May 2005
 
  
% set and check dimensions
nrp = size(dbase,1);
invivosignal = reshape(invivosignal,length(invivosignal),1);
if (nrp~=length(invivosignal))
  warning(['INVIVOSIGNAL has different length than metabolites in ' ...
           'data base. Fixed by trancating / adding zeros.'])
  invivosignal = [invivosignal; zeros(nrp,1)];
  invivosignal = invivosignal(1:nrp);
end
nrbasis = size(dbase,2);

% check arguments and set defaults
if (nargin < 3) | (isempty(t)), t = 1:nrp; end
if (nargin < 4) | (isempty(step)), step = 1; end
if (nargin < 5) | (isempty(abs_scl_factor)), 
  abs_scl_factor = ones(nrbasis,1); end
if (nargin < 6) | (isempty(demp0)),  demp0 = zeros(nrbasis,1); end
if (nargin < 7) | (isempty(freq0)),  freq0 = zeros(nrbasis,1); end
if (nargin < 8) | (isempty(gauss0)), gauss0 = zeros(nrbasis,1); end
if (nargin < 9) | (isempty(e0)),     e0 = zeros(nrbasis,1); end
if (nargin < 10) | (isempty(loline)), loline = 1; end
if (nargin < 11) | (isempty(galine)), galine = 0; end
if (nargin < 12) | (isempty(voline)), voline = 0; end
if (nargin < 13) | (isempty(phasedistort)), phasedistort = 0; end
if (nargin < 14) | (isempty(filtertype)), filtertype = 'none'; end
if (nargin < 15) | (isempty(flow)), flow = -1; end
if (nargin < 16) | (isempty(fhigh)), fhigh = 1; end
if (nargin < 17) | (isempty(ripple)), ripple = 0.01; end
if (nargin < 18) | (isempty(filterlength)), filterlength = 1; end
if (nargin == 19) 
  error('Both PRIOR_EQUAL and EQUAL_TO should be given.')
elseif (nargin < 20) | (length(prior_equal)~=length(equal_to))
  error('PRIOR_EQUAL and EQUAL_TO have different lengths.')
end
if (nargin < 21) | (isempty(allow_change)), allow_change = 0.1; end
if (nargin < 22) | (isempty(baseline_boolean)), baseline_boolean = 0; end
if (nargin < 23) | (isempty(lambda)), lambda = 0.01; end
if (nargin < 24) | (isempty(maxiter)), maxiter = 10; end
if (nargin < 25) | (isempty(plotting)), plotting = 1; end


if (baseline_boolean ~= 1)
  
  error(['When baseline is not required, use AKSES_FIT instead of ' ...
         'AKSES_GCV.']);
  return
end

% initialize baseline-related variables
sord  = 2;                % order of splines
bdeg  = 3;                % degree of splines
li    = 10;               % interval length
ndx   = ceil(nrp/li);     % number of inner knots
nosplpars = ndx + bdeg;    % number of spline pars
xlow  = nrp*(flow+1);
xhi   = nrp*(fhigh+1);
xtemp = min(xlow,xhi);
xhi   = max(xlow,xhi);
xlow  = xtemp;

gcv = []; 
baseline_information = init_baseline(xlow,xhi,nrp,ndx,bdeg,sord,li,lambda);

reg_param = 10.^[0:-.5:-3];
npoints = length(reg_param);

for i = 1:npoints
  
  lambda = reg_param(i);

  if i>1,
  
  % save previously computed information
  [original,recon,baseline,filrecon,filbaseline,ampl,demp,phas,freq,gauss,e,...
   ampl_ratios,ampl_absolute,CR_error,individual_components,residual] = ...
      save_previous(original1,recon1,baseline1,filrecon1,filbaseline1,...
                    ampl1,demp1,phas1,freq1,gauss1,e1,ampl_ratios1,...
                    ampl_absolute1,CR_error1,individual_components1,residual1);
  end
  
  % call akses' fitting function
  [recon1,baseline1,filrecon1,filbaseline1,original1,filter_information,...
   ampl1,demp1,phas1,freq1,gauss1,e1,ampl_ratios1,ampl_absolute1,...
   CR_error1,individual_components1,residual1] = ...
      akses_fit(invivosignal,dbase,t,step,abs_scl_factor,...
                demp0,freq0,gauss0,e0,loline,galine,voline,phasedistort,...
                filtertype,flow,fhigh,ripple,filterlength,...
                prior_equal,equal_to,allow_change,...
                baseline_boolean,lambda,maxiter,plotting);

 
  % compute GCV score
  gcv = [gcv; gcrossval(residual1,baseline_information.alpha,...
                        baseline_information.beta,lambda,...
                        nrp,baseline_information.nosplpars)];

  % test if GCV score start increasing
  if min(gcv)<gcv(end), lambda_opt = reg_param(i-1), break; end
  
  % if GCV didn't reach minimum, set the last lambda as optimal
  if (i == npoints)
    
    lambda_opt = lambda
    
    [original,recon,baseline,filrecon,filbaseline,...
     ampl,demp,phas,freq,gauss,e,ampl_ratios,ampl_absolute,...
     CR_error,individual_components,residual] = ...
        save_previous(original1,recon1,baseline1,filrecon1,filbaseline1,...
                      ampl1,demp1,phas1,freq1,gauss1,e1,ampl_ratios1,...
                      ampl_absolute1,CR_error1,individual_components1,...
                      residual1);

  end
  
end

