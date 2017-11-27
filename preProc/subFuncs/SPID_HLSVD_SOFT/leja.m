function [x_out] = leja(x_in)
%
%    Program orders the values x_in (supposed to be the roots of a 
%    polynomial) in this way that computing the polynomial coefficients
%    by using the m-file poly yields numerically accurate results.
%    Try, e.g., 
%               z=exp(j*(1:100)*2*pi/100);
%    and compute  
%               p1 = poly(z);
%               p2 = poly(leja(z));
%    which both should lead to the polynomial x^100-1. You will be
%    surprised!
%
%    leja is by Markus Lang, and is available from the
%    Rice University DSP webpage: www.dsp.rice.edu




x = x_in(:).'; n = length(x);

a = x(ones(1,n+1),:);
a(1,:) = abs(a(1,:));

[dum1,ind] = max(a(1,1:n));  
if ind~=1
  dum2 = a(:,1);  a(:,1) = a(:,ind);  a(:,ind) = dum2;
end
x_out(1) = a(n,1);
a(2,2:n) = abs(a(2,2:n)-x_out(1));

for l=2:n-1
  [dum1,ind] = max(prod(a(1:l,l:n)));  ind = ind+l-1;
  if l~=ind
    dum2 = a(:,l);  a(:,l) = a(:,ind);  a(:,ind) = dum2;
  end
  x_out(l) = a(n,l);
  a(l+1,(l+1):n) = abs(a(l+1,(l+1):n)-x_out(l));
end
x_out = a(n+1,:);