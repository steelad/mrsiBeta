function [fir_h,h1] = minphJB(h0)
% Transforms linear-phase low-pass filter into a minimum-phase filter

M  = max(size(h0));
rh = roots(h0);
ah = abs(rh);
k  = 1;
k2 = 1;
for t = 1:M-1
    if ((ah(t) > 1))
    rn(k) = rh(t)^-1;
    k = k+1;
    else 
        rn2(k2)=rh(t);
        k2 = k2+1;
  end
end
if exist('rn')
    h1 = poly(leja(rn));
else
    h1 = 1;
end

if exist('rn2')
    h2 = poly(leja(rn2));
else
    h2 = 1;
end
h3 = conv(h1,h2);
[a b] = freqz(h3,1,linspace(0,pi/40,50));
amp = mean(abs(a));
fir_h = h3/amp;
