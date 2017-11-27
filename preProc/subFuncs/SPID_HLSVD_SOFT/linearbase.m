function c = linearbase(data,dbase,ampl,demp,gauss,freq,phas,t,...
                        filter_information,baseline_information)
% linear least-squares fit of linear parameters (amplitudes and
% phases) for given nonlinear parameters (demp,gauss,freq)
  

%filter matrix and data

if ((nargin<8)|(~isempty(baseline_information)))

  j = sqrt(-1);
  [nrp,nrbasis] = size(dbase);

  recon = zeros(nrp,1);
  for i = 1:(nrbasis)
    recon = recon + ampl(i)*dbase(1:nrp,i).*...
            exp((-demp(i)-gauss(i).*t.'+j*freq(i)).*t.'+j*phas(i));
  end

  filrecon = filtering(recon,filter_information,t);

  B = baseline_information.B;
  D = baseline_information.D;
  nosplpars = baseline_information.nosplpars;
  lambda = baseline_information.lambda;
  Btemp = filtering(B,filter_information,t);
    
  Btemp = [Btemp; lambda*D];
  filrecon = [filrecon;zeros(size(baseline_information.D,1),1)];

  c = Btemp\(data-filrecon);
  
else
    c = [];
end


