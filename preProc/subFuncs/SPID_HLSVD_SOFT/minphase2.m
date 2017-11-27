% function y = minphase(x)
% %minphase - minimum phase from mixed phase FIRs
% % y = minphase(x) returns the minimum phase versions of the FIRs
% % in the lines of x
% 
% x = x';
% %Double the impulse response lengths with zeros
% x = [x; zeros(size(x))];
% %Number of impulse responses
% m = size(x,2);
% %Size of each
% n = size(x,1);
% %The weighting function
% wn = [ones(1,m); 2*ones(n/2-1,m); ones(1,m); zeros(n/2-1,m)];
% %The work
% y = real(ifft(exp(fft(wn.*real(ifft(log(abs(fft(x)))))))));
% %Return half the result
% y = y(1:n/2);
% y = y';


function y = minphase(x)
%MINPHASE minimum phase
% MINPHASE(x) returns minumum phase versions of impulse responses in x
% y = MINPHASE(x) returns the minimum phase impulse responses derived from x.
% x = column matrix of impulse responses
x = x';
m = size(x,2);
n = size(x,1);
odd = fix(rem(n,2));
wn = [ones(1,m); 2*ones((n+odd)/2-1,m) ; ones(1-rem(n,2),m);
zeros((n+odd)/2-1,m)];
y = real(ifft(exp(fft(wn.*real(ifft(log(abs(fft(x)))))))));
y = y';