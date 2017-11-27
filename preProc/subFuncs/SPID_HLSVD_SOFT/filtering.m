function filteredsignal = filtering(signal,filter_information,t)
  
% outputs a filtered signal, given an original signal and
% information about the filter to be used
  
if (strcmp(filter_information.type,'fir') == 1) | ...
      (strcmp(filter_information.type,'firmex') == 1) 
    
  fir_h = filter_information.fir_h;
  
  if ~(length(fir_h) == 1 & fir_h(1) == 1)
    filteredsignal = filter(flipud(fir_h(:)),1,signal);
    filteredsignal = filteredsignal(length(fir_h):end,:);
  else
    filteredsignal = signal;
  end
    
elseif (strcmp(filter_information.type,'hlsvd')==1)

  Kest = 50;
  fsampl = 1000;
  time = t/fsampl;
  M = round(length(signal)/2);
  if (size(signal,1)>1) & (size(signal,2)>1)
    filteredsignal = zeros(size(signal)); Kest = 15;
    for i = 1:size(signal,2)
      freqrange = [filter_information.flow filter_information.fhigh];  
      [filteredsignal(:,i),freq,damp,ampl,phas] = HLSVDPROrec(signal(:,i),Kest,fsampl,time,M,freqrange);
    end
  else
    freqrange = [filter_information.flow filter_information.fhigh];  
    filteredsignal = HLSVDPROrec(signal,Kest,fsampl,time,M,freqrange);
  end

elseif  (strcmp(filter_information.type,'fir0') == 1)
    fir_h = filter_information.fir_h;
  
  if ~(length(fir_h) == 1 & fir_h(1) == 1)
    filteredsignal = filter(fir_h(:),1,signal);
    filteredsignal = filteredsignal(length(fir_h):end,:);
  else
    filteredsignal = signal;
  end
    
else
  
  filteredsignal = signal;
  
end
  
  