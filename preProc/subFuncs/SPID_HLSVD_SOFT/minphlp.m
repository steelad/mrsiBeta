function fir_h = minphlp(h0)
% Transforms linear-phase low-pass filter into a minimum-phase filter

M  = max(size(h0));
rh = roots(h0);
ah = abs(rh);
k  = 1;
for t = 1:M-1
  if ((ah(t) < 0.995))
    rn(k) = rh(t);
    k = k+1;
  end
end
k=1;
for t=1:M-1
  if ((ah(t) > 0.99) & (ah(t)<1.01))
    if (ah(t)>1)
      rn2(k) = rh(t)/ah(t)^2;
    else
      rn2(k) = rh(t);
    end
    k=k+1;
  end
end

h1 = poly(rn);
h3 = poly(rn2);
h2 = conv(h1,h1);
h4 = conv(h2,h3);
[a b] = freqz(h4,1,linspace(0,pi/40,50));
amp = mean(abs(a));
fir_h = h4/amp;