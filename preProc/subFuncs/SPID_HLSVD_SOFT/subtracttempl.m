%***************************************************************************
%                                SUBTRACTTEMPL.M
%***************************************************************************
% PURPOSE:  subtract method using AQSES
%***************************************************************************
% CALL:     [scores,misc] = subtract(signal,step,frequency,ndp,begin,trunc,bascorr,dbase,abs_scl_factor,...
%    lineshape, phasedistort,filtertype,boundL,boundH,ripple,filterlength,...
%    prior_equal,equal_to,allow_damp,allow_freq,equalph,baseline_boolean,lambda,...
%    maxiter,plotting,wname,level,adl,thr,sorh,keepapp)
%***************************************************************************
% INPUT:    signal      -- signal                        row vector
%           step        -- time step between points (ms) scalar
%           frequency   -- spectrometer frequency (kHz)  scalar
%           ndp         -- number of data points         scalar
%           begin       -- begin time (ms)               scalar
%           trunc       -- number of truncated points to scalar
%                       estimate the baseline
%           bascorr     -- baseline fitting method       scalar
%                       1:HLSVD-PRO, 2:wavelet (denoising)
%           iter        -- 1:iterations,0: no iter.      binary  
%           dbase       -- database name (met. profiles) string cell
%           abs_scl_factor-- scaling factor              scalar
%           lineshape   -- 1=Lorentzian,2=Gaussian,      scalar
%                           3 = Voigt
%           phasedistort-- phase distortion              binary
%           filtertype  -- 'fir' for MP-FIR              string
%                             (Sundin et al.)
%                          'hlsvd' for HLSVD-PRO
%                             (Laudadio et al.)
%           boundL      -- low bound                     scalar
%           boundH      -- high bound                    scalar
%           ripple      -- ripples (FIR filter)          scalar
%           filterlength-- nbr of filter coefficients	 scalar
%           prior_equal --                               vector
%           equal_to    -- vectors of indeces specifying vector
%                       simple prior knowledge, i.e., dampings
%                       and frequences of metabolites in prior_equal
%                       are equal to equal_to (default [], [])
%           allow_damp  -- allowed damping variations    scalar
%           allow_freq  -- allowed frequency variations	 scalar
%           equalph     -- 1 = equal phase               binary
%                          0 = non equal phase
%           baseline_boolean-- 1 = baseline in the model binary
%                              0 = no basline in the model
%           lambda      -- regularization parameter      scalar
%           maxiter     -- maximum nbr of iterations     scalar
%           plotting    -- display (=1) or not (=0)      binary
%           wname       -- wavelet name                  cell of 1 string
%           level       -- Level at which the wavelet    scalar
%                       decomposition is performed
%           adl         -- automatic def. of defaults    binary
%           thr         -- Global positive threshold     scalar
%           sorh        -- 1:hard, 0:soft thresholding   binary
%           keepapp     -- 1:approximation coefficients  binary
%                       cannot be thresholded, or 0:not  
% OUTPUT:   scores      -- new scores (best variables)   vector
%           misc        -- parameters (ampl., damp., etc)structure
%                       modeled baseline, residual, etc.
%**************************************************************************
function [scores,misc] = subtracttempl(signal,step,frequency,ndp,begin,trunc,bascorr,dbase,abs_scl_factor,...
    lineshape, phasedistort,filtertype,boundL,boundH,ripple,filterlength,...
    prior_equal,equal_to,allow_damp,allow_freq,equalph,baseline_boolean,lambda,...
    maxiter,plotting,wname,level,adl,thr,sorh,keepapp)

torig = [0:step:(ndp-1)*step];
%% TRUNCATION
signalorig = signal;
signal = signal(:,trunc+1:end);

%% FITTING THE METABOLITES (from the truncated signal)
begin = trunc*step;
t   = [0+begin:step:(ndp-1)*step] ;%in ms

if  filtertype==0
    filtertype = '';
elseif filtertype==1
    filtertype = 'fir';
elseif filtertype==2
    filtertype = 'hlsvd';
end

if lineshape==1 %Lorentzian
    loline = 1;
    galine = 0;
    voline = 0;
elseif lineshape==2%Gaussian
    loline = 0;
    galine = 1;
    voline = 0;
elseif lineshape==3%Voigt
    loline = 0;
    galine = 0;
    voline = 1;
end

fsampl = 1/step;
xmin_mouse = ((boundL-4.7)/(fsampl/(frequency/1000)))/1000;%in kHz
xmax_mouse = ((boundH-4.7)/(fsampl/(frequency/1000)))/1000;%in kHz

a = load(char(dbase));
aN = fieldnames(a);
signaltemp = [];
compNames = {};
if isfield(a,'ndp')&isfield(a,'begin')&isfield(a,'step')&isfield(a,'frequency')
else
    error('The variables ndp, begin, step and frequency must exist')
    return
end
k=1;
for i=1:length(aN)
    if (a.ndp == size(a.(aN{i}),2))
        signaltemp = [signaltemp;a.(aN{i})];
        compNames(k) = {aN{i}};
        k = k+1;
    end
end
dbase_landmark1 = signaltemp; % each metabolite profile is a column
%dbaseorig = dbase_landmark1;
sdb = size(dbase_landmark1,2);
if sdb~=ndp
        diff = ndp-sdb;
        if diff<0
            warning('The database has been truncated to fit the number of points of the signal')
            dbase_landmark1 = dbase_landmark1(:,1:ndp);
        else
            warning('The database has been extended by', diff,' to fit the number of points of the signal')
            dbase_landmark1 = [dbase_landmark1 zeros(size(dbase_landmark1,1),diff)];
        end
end
dbaseorigC = dbase_landmark1;
dbase_landmark1 = dbase_landmark1(:,trunc+1:end)'; % each metabolite profile is a column
j=sqrt(-1);
dbase_landmark1 = real(dbase_landmark1)-j*(imag(dbase_landmark1));

%abs_scl_factor = eval(abs_scl_factor);
abs_scl_factor = ones(size(dbase_landmark1,2),1)*abs_scl_factor;
index1 = [1:size(dbase_landmark1,2)];
nrstofs = size(dbase_landmark1,2);
squares = zeros(nrstofs,1);

% initial values for nonlinear parameters
zer    = zeros(size(index1));
demp0  = zer;
gauss0 = zer;
freq0  = zer;
e0     = zer;

%prior knowledge
prior_equal = '';
equal_to = '';

if (lambda<0)
    [recon,baseline,filrecon,filbaseline,filoriginal,filter_information,...
        ampl,demp,phas,freq,gauss,e,ampl_ratios,ampl_absolute,CR_error,...
        individual_components,residual] = ...
        akses_gcv(signal,dbase_landmark1,t,step,abs_scl_factor,...
        demp0,freq0,gauss0,e0,loline,galine,voline,phasedistort,...
        filtertype,xmin_mouse,xmax_mouse,ripple,filterlength,...
        prior_equal,equal_to,allow_damp,allow_freq,baseline_boolean,lambda,...
        maxiter,plotting);
else
    [recon,baseline,filrecon,filbaseline,filoriginal,filter_information,...
        ampl,demp,phas,freq,gauss,e,ampl_ratios,ampl_absolute,...
        CR_error,individual_components,residual] = ...
        akses_fit(signal,dbase_landmark1,t,step,abs_scl_factor,...
        demp0,freq0,gauss0,e0,loline,galine,voline,phasedistort,...
        filtertype,xmin_mouse,xmax_mouse,ripple,filterlength,...
        prior_equal,equal_to,allow_damp,allow_freq,equalph,baseline_boolean,lambda,...
        maxiter,plotting);
end
% recon = reshape(recon,ndp,1);
% baseline = reshape(baseline,ndp,1);
 nn = max(size(filrecon));
% filrecon = reshape(filrecon,nn,1);
% filbaseline = reshape(filbaseline,nn,1);
% filoriginal =  reshape(filoriginal,nn,1);
% individual_components = reshape(individual_components,length(compNames),nn);

scores = ampl;
misc.recon = (real(recon)-sqrt(-1)*imag(recon))'; % reconstructed signal 
misc.baseline =(real(baseline)-sqrt(-1)*imag(baseline))'; %modeled baseline 
misc.filrecon = (real(filrecon)-sqrt(-1)*imag(filrecon))'; %reconstructed filtered signal 
misc.filbaseline = (real(filbaseline)-sqrt(-1)*imag(filbaseline))'; %filtered baseline
misc.filoriginal = (real(filoriginal)-sqrt(-1)*imag(filoriginal))'; %filtered original signal
misc.CR_error = CR_error; % Cramer Rao Bounds
misc.individual_components = individual_components; % corrected metabolite profiles (from the in vitro or simulated database)
misc.residual = residual'; 
misc.ampl = ampl; % amplitude estimates
misc.freq = freq/2/pi; % frequency estimates
misc.damp = demp; % damping estimates 
misc.phas = phas; % phase estimates
misc.begin = begin; % begin point
misc.step = step; % step 
misc.frequency = frequency; % sectrometer frequency
misc.ndp = ndp; %number of data points
misc.compNames = compNames; % name of the metabolites in the database of metabolite profiles

% procres.scores = scores;
% procres.misc = misc;
% save('test.mat','procres')


%% BASELINE FITTING (from the  original signal - reconstructed metabolite
%signal)
j = sqrt(-1);
metabrecon=zeros(size(torig));
for i = 1:nrstofs
  signal_recon = ampl(i)*dbaseorigC(i,:).*(exp((-demp(i)-gauss(i).*torig+j*freq(i)).*torig+j*phas(i)));
  metabrecon = metabrecon + signal_recon;
end

% figure
% plot(real(fftshift(fft(signalorig))))
% hold on;
% plot(real(fftshift(fft(metabrecon))),'r')

basenoise = signalorig-metabrecon; % baseline+noise

% figure
% plot(real(fftshift(fft(basenoise))))

if bascorr==1 %using HLSVD-PRO to model the baseline 
    M = 15;
    MM = round(ndp/2);
    xmin = (boundL-4.7)*frequency/10^6;%in kHz
    xmax = (boundH-4.7)*frequency/10^6;%in kHz
    freqrange1 = [-(fsampl*1000/2-1),xmin*1000]; %in Hz
    freqrange2 = [xmax*1000,(fsampl*1000/2-1)]; %in Hz
    freqrange = [freqrange1;freqrange2];
    [s,freq,damp,ampl,phas,SV] = HLSVDPROrec(basenoise,M,fsampl,torig/1000,MM,freqrange);
    baseRecon=reconsd(torig/1000,freq,damp,ampl,phas); %reconstructed signals    
else %using wavelets to estimate the baseline by denoising the "baseline+noise" signal
    baseRecon = denoisWav(basenoise,step,frequency,ndp,begin,wname,level,adl,thr,sorh,keepapp,plotting);
end

% figure
% plot(real(fftshift(fft(baseRecon))))

pause
%% ESTIMATING THE METABOLITES
metab = signalorig-baseRecon;

% figure
% plot(real(fftshift(fft(signalorig))))
% hold on;
% plot(real(fftshift(fft(metab))),'r')
% plot(real(fftshift(fft(dbaseorigC*misc.ampl))),'g')

dbaseorigC = dbaseorigC'; % each metabolite profile is a column
j=sqrt(-1);
dbaseorigC = real(dbaseorigC)-j*(imag(dbaseorigC));

%abs_scl_factor = eval(abs_scl_factor);
%abs_scl_factor = ones(size(dbase_landmark1,2),1)*abs_scl_factor;
index1 = [1:size(dbase_landmark1,2)];
nrstofs = size(dbase_landmark1,2);
squares = zeros(nrstofs,1);


if (lambda<0)
    [recon,baseline,filrecon,filbaseline,filoriginal,filter_information,...
        ampl,demp,phas,freq,gauss,e,ampl_ratios,ampl_absolute,CR_error,...
        individual_components,residual] = ...
        akses_gcv(metab,dbaseorigC,torig,step,abs_scl_factor,...
        demp0,freq0,gauss0,e0,loline,galine,voline,phasedistort,...
        filtertype,xmin_mouse,xmax_mouse,ripple,filterlength,...
        prior_equal,equal_to,allow_damp,allow_freq,baseline_boolean,lambda,...
        maxiter,plotting);
else
    [recon,baseline,filrecon,filbaseline,filoriginal,filter_information,...
        ampl,demp,phas,freq,gauss,e,ampl_ratios,ampl_absolute,...
        CR_error,individual_components,residual] = ...
        akses_fit(metab,dbaseorigC,torig,step,abs_scl_factor,...
        demp0,freq0,gauss0,e0,loline,galine,voline,phasedistort,...
        filtertype,xmin_mouse,xmax_mouse,ripple,filterlength,...
        prior_equal,equal_to,allow_damp,allow_freq,equalph,baseline_boolean,lambda,...
        maxiter,plotting);
end
signal = recon;


%% RESULTS

scores = ampl;
misc.recon = (real(recon)-sqrt(-1)*imag(recon))'; % reconstructed signal 
%misc.baseline =(real(baseline)-sqrt(-1)*imag(baseline))'; %modeled baseline 
misc.baseline = baseRecon;
misc.filrecon = (real(filrecon)-sqrt(-1)*imag(filrecon))'; %reconstructed filtered signal 
misc.filbaseline = (real(filbaseline)-sqrt(-1)*imag(filbaseline))'; %filtered baseline
misc.filoriginal = (real(filoriginal)-sqrt(-1)*imag(filoriginal))'; %filtered original signal
misc.CR_error = CR_error; % Cramer Rao Bounds
misc.individual_components = individual_components; % corrected metabolite profiles (from the in vitro or simulated database)
misc.residual = residual'; 
misc.ampl = ampl; % amplitude estimates
misc.freq = freq/2/pi; % frequency estimates
misc.damp = demp; % damping estimates 
misc.phas = phas; % phase estimates
misc.begin = begin; % begin point
misc.step = step; % step 
misc.frequency = frequency; % sectrometer frequency
misc.ndp = ndp; %number of data points
misc.compNames = compNames; % name of the metabolites in the database of metabolite profiles

