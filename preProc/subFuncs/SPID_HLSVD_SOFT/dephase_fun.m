function sp_dp = dephase_fun(sp_ft,phc0,phc1)
%function sp_dp = dephase_fun(sp_ft,phc0,phc1)
%written by Chen Li 
%input:   sp_ft - complex data after fourier transform
%         phc0  - the zero order phase correction
%         phc1  - the first order phase correction
%output:  sp_dp - spectral data after phase correction

phc0=phc0*pi/180;              %convert degree to radian
phc1=phc1*pi/180;              %convert degree to radian
[m,n]=size(sp_ft);
a_num=[1:n]./n;
a=phc0.*ones(1,n)+a_num.*phc1; % calculate a(i,j)
re=real(sp_ft);                % get the real part of complex numbers
im=imag(sp_ft);                % get the imaginary part of complex numbers
re_new=re.*cos(a)-im.*sin(a);  % get the new real part by phase correction
im_new=re.*sin(a)+im.*cos(a);  % get the new imaginary part by phase correction

sp_dp=re_new+im_new*i;         % dephased nmr spectral data