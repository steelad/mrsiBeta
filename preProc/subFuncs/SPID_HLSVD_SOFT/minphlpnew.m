function [fir_h,h1] = minphlpnew(h0)
% Transforms linear-phase low-pass filter into a minimum-phase filter
global h1 
global h2
global h3 
global h4

M  = max(size(h0));
rh = roots(h0);
% figure
% plot(rh,'o')
ah = abs(rh);
k  = 1;
for t = 1:M-1
  if ((ah(t) < 0.99))
    rn(k) = rh(t);
    k = k+1;
  end
end
k=1;
for t=1:M-1
  if ((ah(t) > 0.99) & (ah(t)<1.01))
      rn2(k) = rh(t);
      k=k+1;
  end
end

%fvtool(rn)
%fvtool(rn2)

h1 = poly(leja(rn));
h3 = poly(leja(rn2));
h1 = real(h1);
h3 = real(h3);
% h1 = poly(rn);
% h3 = poly(rn2);
h2 = conv(h1,h1);
h4 = conv(h2,h3);

% Adjust for the amplitude
% [a b] = freqz(h4,1,linspace(0,pi/8,50));
% amp = mean(abs(a));
% fir_h = h4/amp;
desired = max(abs(freqz(h0)));
actual = max(abs(freqz(h4)));
fir_h = h4 * (desired/actual);

% figure
% plot(roots(fir_h),'o')
