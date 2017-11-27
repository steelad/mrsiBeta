% Akses one run: fit with all metabolites in dbase.
% Script that calls once the 'akses_fit' function.
% - initializes parameters
% - reads signal
% - reads database
% - filters
% - sets prior knowledge
% - calls fitting function
% - plots results

clear all
close all

% initializations

[load_saved_signal,baseline_boolean,allow_change,phasedistort,...
          plotting,signal_type,bsline_type,noise_level,lambda,maxiter] = ...
    set_hyperparameters; %(0); % put a dummy parameter to skip user interaction

lambda = -1;
nrp = 512;
t0  = 0; tstep = 1; step = 1;
t   = [0:tstep:(nrp-1)*tstep] + t0;
volume = 1;

% construct database
[dbase,abs_scl_factor,name_metabolite] = read_dbase(nrp,filtertype);
%dbase = dbase(:,[1,2,24,25]);abs_scl_factor = abs_scl_factor([1,2,24,25]);
nrstofs = size(dbase,2);

% construct simulation signal

if (load_saved_signal == 1)       %load saved signal 
  
    load('testdat/testnl')
    invivosignal = signaltot(1:nrp);
    
elseif (load_saved_signal == 2)   % simulate sum of basis functions

    [invivosignal,ampl_simul,demp_simul,freq_simul,phas_simul] = ...
        simulate_signal(dbase,nrstofs,nrp,t,abs_scl_factor,signal_type);
    if (baseline_boolean == 1)    % simulate baseline 
      baseline_simul = simulate_baseline(bsline_type,nrp); 
      invivosignal = invivosignal + baseline_simul;
    end
    invivosignal = invivosignal + noise_level*randn(1,nrp);
    
else                              % use your preferred signal here
  
  file2read = 'testdat/philips/te23/LEEN_NORMAAL';
  [signal_save, ndp, step, nfids, volume] = read_philips_file(file2read);
  invivosignal = signal_save(1:nrp); 

end


% compute filter
%xmin_mouse = -0.265019;
%xmax_mouse = -0.051088;
%fir_h = pbfirold(xmin_mouse*2*tstep,xmax_mouse*2*tstep,invivosignal,0.001,60);
filtertype = 'fir';
flow  = -0.265019;
fhigh = -0.051088;
ripple = 0.01;
filterlength = 80;

% prior knowledge for equal dampings
[prior_equal,equal_to] = set_prior({'gln','gaba','pch','naag','lip2'},...
                                   {'glu','cr','gpc','naa','lip1'},...
                                   name_metabolite);

% initial values for nonlinear parameters
demp0      = zeros(nrstofs,1);
gauss0     = zeros(nrstofs,1);
freq0      = zeros(nrstofs,1);
e0         = zeros(nrstofs,1);

% call fakses function
  
[recon,baseline,filrecon,filbaseline,original,filter_information,...
 ampl,demp,phas,freq,gauss,e,ampl_ratios,ampl_absolute,CR_error_all,...
 individual_components,residual] = ...
  akses_gcv(invivosignal,dbase,t,step,abs_scl_factor,...
            demp0,freq0,gauss0,e0,1,0,0,phasedistort,...
            filtertype,flow,fhigh,ripple,filterlength,...
            prior_equal,equal_to,allow_change,baseline_boolean,lambda,...
            maxiter,plotting);

CR_error = CR_error_all.a;
% plot fitting results
disp_concentrations

