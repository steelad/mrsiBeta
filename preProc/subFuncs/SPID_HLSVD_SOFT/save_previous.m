function [original,recon,baseline,filrecon,filbaseline,...
          ampl,demp,phas,freq,gauss,e,ampl_ratios,ampl_absolute,...
          CR_error,individual_components,residual] = ...
      save_previous(original0,recon0,baseline0,filrecon0,filbaseline0,...
                    ampl0,demp0,phas0,freq0,gauss0,e0,ampl_ratios0,...
                    ampl_absolute0,CR_error0,individual_components0,residual0);

  
% copy variable0 to variable
original = original0;
recon    = recon0;
baseline = baseline0;
filrecon = filrecon0;
filbaseline = filbaseline0;
ampl     = ampl0;
demp     = demp0;
phas     = phas0;
freq     = freq0;
gauss    = gauss0;
e        = e0;
ampl_ratios = ampl_ratios0;
ampl_absolute = ampl_absolute0;
CR_error = CR_error0;
individual_components = individual_components0;
residual = residual0;