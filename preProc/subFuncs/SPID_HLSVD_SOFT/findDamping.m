% find the damping of the peaks in a specific frequency region of a signal
function damping = findDamping(signal,step,ndp)

%the spectra are supposed to be in absorption mode (maximum real part,
%minimum imaginary part (peaks in the real part are pointing upwards)
msF = fftshift(fft(signal));
[mx,indM] = max(abs(msF));

% ldiff=[msf(1) msf(2:end)-msf(1:end-1)];%left diff
% rdiff=[msf(1:end-1)-msf(2:end) msf(end)];%right diff

sa=round(ndp/20); %searching around the maximum
if indM-sa<1
    e1 = 1;
else
    e1 = indM-sa;
end
if indM+sa>length(signal)
    e2 = length(signal);
else
    e2 = indM+sa;
end
[mn,indm] = min(abs(real(msF(e1:e2))-real(mx)/2));

HH=abs(msF(e1+indm-1)); %half height

xlow  = -1000 / (2*step) ;
xstep = 1000 / (ndp*step);
xupp = -xlow;
xpts     = xlow:xstep:xupp;
kHz = xpts/1000; % in kHz

FD = kHz(indM)-kHz(e1+indm-1); % frequency difference between the highest point of the highest peak and the point at half height of this maximum (in kHz)
damping = abs(2*FD*pi)/(1+exp(mx/2/HH)-exp(1));
if damping<0
    damping=0;
end