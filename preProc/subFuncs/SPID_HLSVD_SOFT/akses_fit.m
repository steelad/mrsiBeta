function [recon,baseline,filrecon,filbaseline,filoriginal,...
          filter_information,ampl,demp,phas,freq,gauss,e,...
          ampl_ratios,ampl_absolute,CR_error,individual_components,...
          residual] = ...
    akses_fit(invivosignal,dbase,t,step,abs_scl_factor,...
              demp0,freq0,gauss0,e0,loline,galine,voline, ...
              phasedistort,filtertype,flow,fhigh,ripple,filterlength,...
              prior_equal,equal_to,allow_change_damp,allow_change_freq,...
              equalph,baseline_boolean,lambda,maxiter,plotting)

% Function that fits a short echo time proton spectrum in the time domain 
% using a metabolite database.
%
% Method: nonlinear least squares with prior knowledge (linear
%         relations between some parameters) and linear bounds.
%
% Input variables:
%   invivosignal  = signal to be quantified
%   dbase         = data base of metabolite signals, stored in columns
%   t             = vector of time instances (default [1:length(invivosignal)])
%   step          = frequency sampling step (default 1)
%   abs_scl_factor= scaling factors for amplitudes (default = ones)
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
%   allow_change_damp = constant value that fixes the linear bounds on
%                       dampings (default 0.01)
%   allow_change_freq = constant value that fixes the linear bounds on
%                       frequencies (default 0.01)
%   equalph           = flag: 0 - free phase corrections (default) 
%                             1 - equal phase corrections for all metab.
%   baseline_boolean  = flag: 0 - no baseline is fit (default) 
%                             1 - baseline is fit
%   lambda = regularization parameter for baseline coefficients
%            (default 0.1)
%   maxiter = maximum number of Levenberg-Marquardt iterations
%             (default 10)
%   plotting = flag: 0 - no plots during optimization process
%                    1 - plot frequency domain fit at every iteration
%
% Output variables:
%   recon    = reconstructed (fitted) signal
%              length:  length(invivosignal);
%   baseline = fitted baseline (time-domain)
%              length:  length(invivosignal);
%   filrecon = filtered reconstructed (fitted) signal
%              length:  length(invivosignal) - length(fir_h) + 1;
%   filbaseline = filtered fitted baseline (time-domain)
%              length:  length(invivosignal) - length(fir_h) + 1;
%   filoriginal = filtered original in vivo signal 
%              length:  length(invivosignal) - length(fir_h) + 1;
%   filter_information = structure containing filter information
%                        (see init_filter.m)
%   ampl, demp, phas, freq, gauss, e = parameters of fitted signal
%              length:  nrbasis = size(dbase,2);
%   ampl_ratios   = amplitude ratios to Cr amplitude
%   ampl_absolute = metabolite concentrations
%   CR_error      = Cramer-Rao error bounds 
%                   a structure containing the fields: a, d, g, f, ph, e,
%                   such that, e.g., CR_error.a is a vector of CR bounds
%                   for the amplitudes
%   individual_components = filtered individual components in data base
%   residual      = residual of the filtered fit

% Diana Sima, KUL-ESAT, May 2005. Modified: Feb. 2007.

% set and check dimensions
nrp = size(dbase,1);
invivosignal = reshape(invivosignal,length(invivosignal),1);
if (nrp~=length(invivosignal))
  warning(['INVIVOSIGNAL has different length than metabolites in ' ...
           'data base. Fixed by truncating / adding zeros.'])
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
if (nargin < 21) | (isempty(allow_change_damp)), allow_change_damp = 0.1; end
if (nargin < 22) | (isempty(allow_change_freq)), allow_change_freq = 0.1; end
if (nargin < 23) | (isempty(equalph)), equalph = 0; end
if (nargin < 24) | (isempty(baseline_boolean)), baseline_boolean = 0; end
if (nargin < 25) | (isempty(lambda)), lambda = 0.01; end
if (nargin < 26) | (isempty(maxiter)), maxiter = 10; end
if (nargin < 27) | (isempty(plotting)), plotting = 1; end

% initializations
demp      = demp0;
gauss     = gauss0;
freq      = freq0;
e         = e0;

% prior knowledge: 
[demptype,gausstype,freqtype,etype] =...
    priorknowledge(nrbasis,prior_equal,equal_to);
                                        
% filter the signal.
filter_information = init_filter(filtertype,flow,fhigh,invivosignal,...
                                 ripple,filterlength,step,size(invivosignal,1));
%filterlength = length(filter_information.fir_h)
%filter_information.fir_h
data = filtering(invivosignal,filter_information,t); 
nrpf = length(data);

% added by JB     
% figure
% disp_sig2(data(:,1).',0,1,0,length(data(:,1)),'',0,0,'nmr',1,0,'r'); %real
% disp_sig2(data(:,1).',0,1,0,length(data(:,1)),'',0,1,'nmr',1,0,'r'); %absolute


% initialize variables to be used in the optimization
if (equalph)
  [ampl,phas] = linear2(data,dbase,demp,gauss,freq,t,filter_information);
else
  ampl = []; phas = [];
end

[alf,ampl,demp,gauss,freq,phas,e,lb,ub,nra,nrd,nrg,nrf,nrph,nre,...
 amplref,dempref,gaussref,freqref,phasref,eref]=...
    init_one_run(ampl,demp,gauss,freq,phas,e,...
                 demptype,gausstype,freqtype,etype,...
                 loline,galine,voline,phasedistort,equalph,...
                 allow_change_damp,allow_change_freq);

   
                 
if (baseline_boolean == 1)   

  % initialize baseline-related variables
  sord  = 2;                 % order of splines
  bdeg  = 3;                 % degree of splines
  li    = 10;                % interval length
  ndx   = ceil(nrp/li);      % number of inner knots
  nosplpars = ndx + bdeg;    % number of spline pars
  xlow  = nrp*(flow+1);
  xhi   = nrp*(fhigh+1);
  xtemp = min(xlow,xhi);
  xhi   = max(xlow,xhi);
  xlow  = xtemp;
  
  % store (in baseline_information) all values to be passed to the solver
  baseline_information = init_baseline(xlow,xhi,nrp,ndx,bdeg,sord,li,lambda);

  % extend data with zeros
  data = [data; zeros(size(baseline_information.D,1),1)];
  if (~equalph)
    dataz = zeros(2*(nrpf-nrbasis-nosplpars+size(baseline_information.D,1)),1);
  else
    dataz = zeros(2*(nrpf-nosplpars+size(baseline_information.D,1)),1);
  end
  
%   whos dataz
%   nrbasis
%   nosplpars
%   ttt = size(baseline_information.D,1)

else 

  % no baseline will be included in the fit
  baseline_information = []; 
  if (~equalph)
    dataz = zeros(2*(nrpf-nrbasis),1);
  else
    dataz = zeros(2*nrpf,1);
  end

end

if (plotting ~= 1), displaying = 'off'; else displaying = 'iter'; end 

% optimization settings
options = optimset('diagnostics','off','Jacobian','on','display',...
                   displaying,'LargeScale','on','LevenbergMarquardt','on',...
                   'MaxIter', maxiter, 'TolFun',1e-6, 'TolX', 1e-6,...
                   'DerivativeCheck','off');

% call the nonlinear least squares solver
disp('NLLS FILTERING')
%whos dbase
[x,resnorm,residual,some_exitflag,some_output,some_lambda,J]  = ...
    lsqcurvefit('funjaco',alf,1,dataz,lb,ub,options,dbase,t,...
                demptype,gausstype,freqtype,etype,amplref,...
                dempref,gaussref,freqref,phasref,eref,filter_information,...
                demp,gauss,freq,e,loline,galine,voline,t(1),...
                phasedistort,equalph,data,filter_information,...
                baseline_information,step,plotting);

% whos residual
            
% extract parameter information from optimal solution x
[ampl,demp,gauss,freq,phas,e] = ...
    update(x,nrbasis,demptype,gausstype,freqtype,etype,...
           amplref,dempref,gaussref,freqref,phasref,eref,demp,gauss,freq,e,...
           loline,galine,voline,phasedistort,equalph);

% compute optimal linear parameters
%disp('LINEAR FILTERING')
if (~equalph)
  [ampl,phas,zlin] = linear2(data,dbase,demp,gauss,freq,t,...
                            filter_information,baseline_information);
else
  zlin = linearbase(data,dbase,ampl,demp,gauss,freq,phas,t,...
                    filter_information,baseline_information);
end

 % prepare plotting
disp('FINAL RESULT FILTERING')
[ampl_ratios,ampl_absolute,recon,filoriginal,filrecon,individual_components,...
 xpts2,hr] = finalresult(ampl,abs_scl_factor,demp,gauss,freq,phas,t,0, ...
                         invivosignal,dbase,filter_information,1,step);
% disp_sig2(recon.',0,1,t(1),length(recon),'',0,0,'nmr',1,0,'g');  % real                                     
filteredrecon = filtering(recon,filter_information,t);
%whos filteredrecon
%disp_sig2(filteredrecon.',0,1,t(1),length(filteredrecon),'',0,1,'nmr',1,0,'g');   %absolute
%load baselinesignal
%disp_sig2(0.05*baseline,0,1,t(1),length(baseline),'',0,1,'nmr',1,0,'k');   %absolute 

if (baseline_boolean == 1)
    zlin = zlin(:);
  baseline = baseline_information.B*zlin;
  baseline = reshape(baseline,size(invivosignal));
  disp('BASELINE FILTERING')
  filbaseline = filtering(baseline,filter_information,t);
  
  
else
  baseline = zeros(size(invivosignal));
  filbaseline = zeros(size(filoriginal));
end

if nargout >= 15
  [CR_error.a,CR_error.d,CR_error.g,CR_error.f,CR_error.ph,CR_error.e] = ...
      estimate_CR_bounds(ampl,demp,freq,phas,gauss,e,t,dbase,...
                         filter_information,filoriginal,filrecon,...
                         filbaseline,loline,galine,voline);
end

