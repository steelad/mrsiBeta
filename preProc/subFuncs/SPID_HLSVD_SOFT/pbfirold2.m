function fir_h = pbfirold(wl,wh,signal,ri,M0)
% Minimum-phase passband FIR filter design.
%       fir_h = pbfir(wl,wh,signal,ri,M0)
% fir_h - FIR filter coeffienct vector
% wl - lower cut-off frequency, Normalized frequency in [-1,1]
% wh - upper cut-off frequency, Normalized frequency in [-1,1]
% signal - complex MRS signal
% ri - passband ripple (typical 0.01)
% M0 - initial filter length (typical 30) 
% Note that filtering in matlab should be done with maximum-phase filter
% and discarding first M-1 samples. Ex.:
% filsignal0=filter(fir_h(end:-1:1),1,signal);
% filsignal=filsignal0(M:end);

N  = max(size(signal));
wc = (wh-wl)/2;

noise = std(real(signal(end-20:end)));
c=1;
while ~noise 
    noise  = std(real(signal(end-20*c:end)));
end
maxs  = max(abs(fft(signal)))/sqrt(N);
sup   = noise/maxs*10;
mnew  = 1e10;
ok = 1;
M  = M0;
%disp('5')
while ok == 1;
%   M 
%   wc
%   ri 
%   sup
  fir_h = fircls1(M,wc,ri,sup);
%   figure
%   plot(fir_h)
  fir_h = minphlp(fir_h);
  M2    = max(size(fir_h));
  fir_h = fir_h.*exp(-j*pi*(wl+wc)*[0:M2-1]);
  f     = filter(fir_h,1,signal(end:-1:1));
  ff    = abs(fftshift(fft(f(end:-1:M),2048)));
  mold  = max(ff(1:round((wl+1)*1024)-10))/sqrt(N);
  mold2 = max(ff(round((wh+1)*1024)+10:end))/sqrt(N);
  mold  = max(mold,mold2);
  if mold < 2*noise;
    ok = 0;
  else
    ok = 1;
    if mold > 0.9*mnew
      M = M+10;
    if M > 80
	ok  = 0;
	sup = supold;
	M   = Mold;
	fir_h = fircls1(M,wc,ri,sup);
	fir_h = minphlp(fir_h);
	M2    = max(size(fir_h));
	fir_h = fir_h.*exp(-j*pi*(wl+wc)*[0:M2-1]);
	fil   = filter(fir_h,1,signal(end:-1:1));
        ff    = abs(fftshift(fft(fil(end:-1:M2-1),2048)));
        mold  = max(ff(1:round((wl+1)*1024)-10))/sqrt(N);
        mold2 = max(ff(round((wh+1)*1024)+10:end))/sqrt(N);
        mold  = max(mold,mold2);
        %disp('end')
      end
      sup = sup*2;
       %disp('6')
    else
      mnew = mold;
      Mold = M;
      supold = sup;
      M = M0;
      end
    sup = sup/2;
    end
end
 