%***************************************************************************
%                                AQSESMRSI.M
%***************************************************************************
% PURPOSE:  Quantitation MRSI data using AQSES_MRSI
%***************************************************************************
% CALL:     [scores,misc] = AQSESMRSI(signal,step,frequency,ndp,begin,position,procres,classopt,dbase,abs_scl_factor,...
%    lineshape,
%    phasedistort,filtertype,boundL,boundH,ripple,filterlength,...
%    prior_equal,equal_to,allow_damp,allow_freq,equalph,baseline_boolean,lambda,...
%    maxiter,plotting,option)
%***************************************************************************
% INPUT:    signal      -- signal                        nos x ndp matrix
%           step        -- time step between points (ms) scalar
%           frequency   -- spectrometer frequency (kHz)  scalar
%           ndp         -- number of data points         scalar
%           begin       -- begin time (ms)               scalar
%           position    -- position of the voxel         matrix
%            matrix of size Lx3 with L= nbr of signals
%            1st column for the row position, 2nd column for the column
%            position of the voxel, 3rd column for the identification of
%            the voxel (from 1 to L).
%           classopt    -- classification options        structure
%               classtype  -- class type                 vector/cell
%           procres     -- processing options            structure
%               scores     -- data (MxN with M=nbr of
%                           data and N=nbr of variables) matrix
%           dbase       -- database name (met. profiles) string cell
%           tissuetype  -- tissue type                   vector of integers
%              e.g.: beta = [1 0 2 1 1 1 0 0 2 1 ... 0]^T  (T for transposed) . beta, here, is just a class labeling vector, and is not
%               similar to beta_cs in the attached document where 0: normal
%               tissue 1: CSF, 2: tumor tissue, for instance
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
%           option      -- quantitation method           scalar
%            0:default just prior knowledge, 1: using spatial info during
%            optimization
%           plotting    -- display (=1) or not (=0)      binary
% OUTPUT:   scores      -- new scores (best variables)   vector
%           misc        -- parameters (ampl., damp., etc)structure
%                       modeled baseline, residual, etc.
%**************************************************************************
function [scores,misc] = AQSESMRSI(signal,step,frequency,ndp,begin,position,procres,classopt,dbase,abs_scl_factor,...
    lineshape, phasedistort,filtertype,boundL,boundH,ripple,filterlength,...
    prior_equal,equal_to,allow_damp,allow_freq,equalph,baseline_boolean,lambda,...
    maxiter,plotting)

'intru aqsesMRSI'
t   = [0:step:(ndp-1)*step] + begin;%in ms

%initialization
if nargin<9
    error('The number of function parameters should be at least equal to 7')
end

if nargin<10
    abs_scl_factor = 1;
end
if nargin<11
    lineshape = 1;
end
if nargin<12
    phasedistort = 0;
end
if nargin<13
    filtertype = 0;
end
if nargin<14
    boundL =0.25;
end
if nargin<15
    boundH =4.2;
end
if nargin<16
    ripple =0.01;
end
if nargin<17
    filterlength =50;
end
if nargin<18
    prior_equal =0;
end
if nargin<19
    equal_to =0;
end
if nargin<20
    allow_damp =0.03;
end
if nargin<21
    allow_freq =0.1885; %in kHz (2*pi*freq)
end
if nargin<22
    equalph =1;
end
if nargin<23
    baseline_boolean=0;
end
if nargin<24
    lambda=1;
end
if nargin<25
    maxiter=25;
end
if nargin<26
    plotting=0;
end
%if nargin<27
 %   option=0;
%end

if exist('classopt')
    classtype=classopt.classtype;
else
    classtype = ones(size(signal,1),1); %we consider by default that all signals have been acquired from tissue of the same type
end

weightfordamp = 0.2;
weightforgauss = sqrt(0.2);
weightforfreq = 2;


if  filtertype==0
    filtertype = '';
elseif filtertype==1
    filtertype = 'fir';
elseif filtertype==2
    filtertype = 'hlsvd';
elseif filtertype==3
    filtertype = 'fir0';
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
sdb = size(dbase_landmark1,2);
if sdb~=ndp
    diff = ndp-sdb;
    if diff<0
        warning('The database has been truncated to fit the number of points of the signal')
        dbase_landmark1 = dbase_landmark1(:,1:ndp);
    else
        warning(['The database has been extended by', diff,' to fit the number of points of the signal'])
        dbase_landmark1 = [dbase_landmark1 zeros(size(dbase_landmark1,1),diff)];
    end
end
dbase_landmark1 = dbase_landmark1'; % each metabolite profile is a column
j=sqrt(-1);
dbase_landmark1 = real(dbase_landmark1)-j*(imag(dbase_landmark1));
% whos dbase_landmark1
%  figure(10)
%  [offset2]=disp_sig2(signal,0,step,begin,length(signal),'',0,0,frequency,1,0);
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
nos = size(signal,1);
opt=0;
beta = classtype2boolean(position,classtype,opt);

if loline
    param = zeros(nos,2*nrstofs);
elseif galine
    param = zeros(nos,2*nrstofs);
else
    param = zeros(nos,3*nrstofs);
end

%% first round, individual voxel parameter estimatian using AQSES

conv=0; % convergence flag (0:no convergence, 1:convergence)
        amplsave = [];
        reconsave = [];
        for i = 1:nos
            [recon,baseline,filrecon,filbaseline,filoriginal,filter_information,...
                ampl,demp,phas,freq,gauss,e,ampl_ratios,ampl_absolute,...
                CR_error,individual_components,residual] = ...
                akses_fit(signal(i,:),dbase_landmark1,t,step,abs_scl_factor,...
                demp0,freq0,gauss0,e0,loline,galine,voline,phasedistort,...
                filtertype,xmin_mouse,xmax_mouse,ripple,filterlength,...
                prior_equal,equal_to,allow_damp,allow_freq,equalph,baseline_boolean,lambda,...
                maxiter,plotting); 
                
            amplsave = [amplsave; ampl];                  
            CR_error_aqses(i) = CR_error;
    if loline
                param(i,:) = [demp freq];
            elseif galine
                param(i,:) = [gauss freq];
            else
                param(i,:) = [demp gauss freq];
    end
            
    misc1(i).recon = (real(recon)-sqrt(-1)*imag(recon))'; % reconstructed signal
    misc1(i).baseline =(real(baseline)-sqrt(-1)*imag(baseline))'; %modeled baseline
    misc1(i).filrecon = (real(filrecon)-sqrt(-1)*imag(filrecon))'; %reconstructed filtered signal
    misc1(i).filbaseline = (real(filbaseline)-sqrt(-1)*imag(filbaseline))'; %filtered baseline
    misc1(i).filoriginal = (real(filoriginal)-sqrt(-1)*imag(filoriginal))';
    misc1(i).ampl = ampl; % amplitude estimates
    misc1(i).freq = freq/2/pi; % frequency estimates
    misc1(i).damp = demp; % damping estimates
    misc1(i).phas = phas*180/pi; % phase estimates in degrees
        end
        save 'amplafter1aqses.mat' amplsave param
        save first_round
       
 
 
%% AQSES-MRSI 
load first_round 
conv=0;

        % Sweeping
        cc=0;
        while conv==0
            paramold = param;
            for i = 1:nos
                [demp0,freq0,gauss0,e0] = medianparam(param(logical(beta(i,:)),:),loline,galine,voline,nrstofs);
                if loline
                    dampmean = mean(param(logical(beta(i,:)),1:nrstofs),1);
                    freqmean = mean(param(logical(beta(i,:)),nrstofs+1:2*nrstofs),1);
                elseif galine
                    freqmean = mean(param(logical(beta(i,:)),nrstofs+1:2*nrstofs),1);
                else
                    dampmean = mean(param(logical(beta(i,:)),1:nrstofs),1);
                    freqmean = mean(param(logical(beta(i,:)),2*nrstofs+1:2*nrstofs),1);
                end

                allow_damp_new = -max(dampmean-.25/(cc+1)*allow_damp,0);
                allow_damp_new = [allow_damp_new,min(dampmean+.25/(cc+1)*allow_damp,allow_damp)];
                allow_freq_new = -max(freqmean-.25/(cc+1)*allow_freq,-allow_freq);
                allow_freq_new = [allow_freq_new,min(freqmean+.25/(cc+1)*allow_freq,allow_freq)];
                
                sigma = std(signal(i,end-20:end))/sqrt(2); %std of the noise in the time domain
            
                
                [recon,baseline,filrecon,filbaseline,filoriginal,...
                    filter_information,ampl,demp,phas,freq,gauss,e,...
                    ampl_ratios,ampl_absolute,CR_error,individual_components,...
                    residual] = ...
                    akses_fitMRSI(signal(i,:),dbase_landmark1,t,step,abs_scl_factor,...
                    demp0,freq0,gauss0,e0,loline,galine,voline, ...
                    phasedistort,filtertype,xmin_mouse,xmax_mouse,ripple,filterlength,...
                    prior_equal,equal_to,allow_damp_new,allow_freq_new,...
                    equalph,baseline_boolean,lambda,maxiter,plotting,...
                    param,beta(i,:),weightfordamp,weightforgauss,weightforfreq,sigma);
                recon = reshape(recon,ndp,1);
                baseline = reshape(baseline,ndp,1);
                nn = max(size(filrecon));
                filrecon = reshape(filrecon,nn,1);
                filbaseline = reshape(filbaseline,nn,1);
                filoriginal =  reshape(filoriginal,nn,1);
                individual_components = reshape(individual_components,length(compNames),nn);

                scores(i,:) = ampl;
                               
                misc(i).recon = (real(recon)-sqrt(-1)*imag(recon))'; % reconstructed signal
                misc(i).baseline =(real(baseline)-sqrt(-1)*imag(baseline))'; %modeled baseline
                misc(i).filrecon = (real(filrecon)-sqrt(-1)*imag(filrecon))'; %reconstructed filtered signal
                misc(i).filbaseline = (real(filbaseline)-sqrt(-1)*imag(filbaseline))'; %filtered baseline
                misc(i).filoriginal = (real(filoriginal)-sqrt(-1)*imag(filoriginal))'; %filtered original signal
                misc(i).CR_error = CR_error; % Cramer Rao Bounds
                misc(i).individual_components = individual_components; % corrected metabolite profiles (from the in vitro or simulated database)
                misc(i).residual = residual';
                misc(i).ampl = ampl; % amplitude estimates
                misc(i).freq = freq/2/pi; % frequency estimates
                misc(i).damp = demp; % damping estimates
                misc(i).phas = phas*180/pi; % phase estimates in degrees
                misc(i).begin = begin; % begin point
                misc(i).step = step; % step
                misc(i).frequency = frequency; % sectrometer frequency
                misc(i).ndp = ndp; %number of data points
                misc(i).compNames = compNames; % name of the metabolites in the database of metabolite profiles
                CR_error_aqses5(i) = CR_error;
                
                if loline
                    param(i,:) = [demp freq];
                elseif galine
                    param(i,:) = [gauss freq];
                else
                    param(i,:) = [demp gauss freq];
                end
            end
            C = 1/(size(param,1)*size(param,2))*sum(sum((paramold-param).^2./param.^2)); % convergence criterion
                                
            if C<1E-3 | cc>10
                conv=1;
            end
            cc = cc+1;           
        end
        save amplafterAqsesMRSI scores CR_error_aqses5 misc 
%

 

