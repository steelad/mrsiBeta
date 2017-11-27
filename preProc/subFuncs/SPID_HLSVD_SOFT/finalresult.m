function [ampl_ratios,ampl_absolute,recon,filoriginal,filrecon,...
          individual_components,xpts2,hr] = ...
    finalresult(ampl,abs_scl_factor,demp,gauss,freq,phas,t,phi0,data,...
               stof,filter_information,volume,step)

% Function that computes reconstructed signal and returns variables
% useful for plotting. 
%
% Output vars:
% 
% ampl_ratios = amplitudes scaled by Cr, using absolute scaling factors
% ampl_absolute = concentrations
% recon = reconstructed signal
% filoriginal = original signal (filtered)
% filrecon = reconstructed signal (filtered)
% individual_components = database components, corrected and filtered
% signaltot = reconstructed signal (unfiltered)
% xpts2 = abscissas used for plotting
% hr = filter used for plotting


j = sqrt(-1);
[nrp,nrstofs] = size(stof);
% initialize total signal to be reconstructed
recon = zeros(nrp,1);
%whos recon
for i = 1:nrstofs
  signal_recon = ampl(i)*stof(:,i).*(exp((-demp(i)-gauss(i).*t+j* ...
                                          freq(i)).*t+j*phas(i))).';
  fildata      = filtering(signal_recon,filter_information,t);
  recon        = recon + signal_recon;
  individual_components (i,:) = fildata(:).';
end

filrecon    = filtering(recon,filter_information,t);

filoriginal = filtering(data,filter_information,t);

% calculate ../Cr ratios based on scaling factor in dbase,
temp_ampl=(ampl(:).*abs_scl_factor(:));
ampl_ratios=temp_ampl/temp_ampl(1);

% rescale absolute values to Creatine==8mM to get rid of possible 
% transmit ampl differences
ampl_absolute = temp_ampl*50*8/volume; % volume voxels in dbase is
                                       % 20x20x20mm= ..= 8 ml. 
				       % note that the value of "50"
                                       % pops up and is hard-coded,
                                       % simply because 50mM is the
                                       % concentration of Creatine
                                       % in the phantom measurements

% variables necessary for plotting
ldata = length(filoriginal);
[xpts,xpts2]  = absi(ldata,step);
filter_information.fir_h =[];
if (isfield(filter_information,'fir_h'))
  fir_h = filter_information.fir_h;
else 
  fir_h = 1;
end
hr =1;
% [hr,w]        = freqz(fir_h,1,ldata,'whole');
% hr(ldata+1)   = hr(1);

