function filter_information = init_filter(type,flow,fhigh,signal,ripple,M,step,ndp)

% initializes information about the filter that should be used
  
filter_information.type = type;

if (strcmp(type,'fir')==1)
  
  filter_information.fir_h = pbfirnew(2*flow,2*fhigh,signal,ripple,M);
  
elseif (strcmp(type,'firmex')==1)

  filter_information.fir_h = pbfirmex(2*fhigh,2*flow,signal,ripple,M);
  filter_information.fir_h = filter_information.fir_h(end:-1:1);
    
elseif (strcmp(type,'hlsvd')==1)

  filter_information.flow  = flow*1000;
  filter_information.fhigh = fhigh*1000;
  
elseif (strcmp(type,'fir0')==1)
  PB = [flow,fhigh];  
  filter_information.fir_h = FIRwat0(PB,signal,M,step,ndp); 
else
  
  % no filter

end
