function res = times(a,b)

% b : susceptibility map in frequency domain

if isa(a,'FT_v2') == 0
    error('In  A.*B only A can be Fourier operator');
end

if a.adjoint
    res = sqrt(length(b(:)))*fftshift(ifftn(ifftshift(b)));
else
    res = 1/sqrt(length(b(:)))*fftshift(fftn(ifftshift(b))) .* a.mask;
end

end

