function plotSpecFinish(inputArray,Vx,Vy,FFT)
if FFT == 1
    temp = squeeze(real(fftshift(fft(conj(inputArray(Vx,Vy,:))))));
else
    temp = squeeze(real(inputArray(Vx,Vy,:)));
end
plot(temp)
xlim([512 1024])
ylim([min(temp) max(temp)]) 