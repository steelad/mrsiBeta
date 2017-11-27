function [demp0,freq0,gauss0,e0] = medianparam(parami,loline,galine,voline,nrstofs)

demp0 = []; freq0 = []; gauss0 = []; 
medianpar = median(parami,1); 
if loline
    demp0 = medianpar(1:nrstofs);
    freq0 = medianpar(nrstofs+1:2*nrstofs);
elseif galine
    gauss0 = medianpar(1:nrstofs);
    freq0  = medianpar(nrstofs+1:2*nrstofs);
else
    demp0 = medianpar(1:nrstofs);
    gauss0 = medianpar(nrstofs+1:2*nrstofs);
    freq0 = medianpar(2*nrstofs+1:2*nrstofs);
end
e0 = zeros(1,nrstofs);