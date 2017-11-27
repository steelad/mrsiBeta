function [ampl,phas,c] = linear2(data,dbase,demp,gauss,freq,t,...
    filter_information,baseline_information)
% linear least-squares fit of linear parameters (amplitudes and
% phases) for given nonlinear parameters (demp,gauss,freq)

j = sqrt(-1);
[nrp,nrbasis] = size(dbase);
for i = 1:(nrbasis)
    A(1:nrp,i) = dbase(1:nrp,i).*exp((-demp(i)-gauss(i).*t.'+j*freq(i)).*t.');
end


%filter matrix and data

Atemp = filtering(A,filter_information,t);

if ((nargin>=8)&(~isempty(baseline_information)))

    B = baseline_information.B;
    D  = baseline_information.D;
    nosplpars = baseline_information.nosplpars;
    lambda = baseline_information.lambda;
    Btemp = filtering(B,filter_information,t);

    Atemp = [Atemp Btemp; zeros(size(D,1),nrbasis) lambda*D];
end

data = reshape(data,size(Atemp,1),1);
c = Atemp\data;

ampl = abs(c(1:nrbasis));
phas = angle(c(1:nrbasis));
c = c(nrbasis+1:end);
%JB correction (s.t. ampl, phas idem as in linearbase)
ampl = ampl';
phas = phas';
c = c';

