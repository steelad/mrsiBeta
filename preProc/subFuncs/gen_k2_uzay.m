function [k]= gen_k2_uzay(n_rings, diff_rings, ntheta, FOV, alpha)

%if nargin < 4
    alpha   =   1;
    beta=1.;
%end
kmax1 = alpha*(n_rings*2)/(2*FOV);
dx= FOV/(n_rings*2);
kmax    =   alpha/(2*dx);
tmp     =   linspace(0,kmax,32768);
n_rings=n_rings*diff_rings;
y       =   (dx*beta/2)*(1+cos(2*pi*tmp*dx/alpha));
y       =   (n_rings)*(y.*tmp)/sum(y.*tmp);

r(1)    =   0;
for i = 1:n_rings-1
   r(i+1)    =   tmp(find(cumsum(y)>=i,1));
end
r(n_rings+1)  =   kmax;
dr  =   diff(r);
r   =   r(1:end-1)+(2/3)*dr;

%r   =   r*pi*(2*dx);

dtheta  =   2*pi/ntheta;

k   =   [];
for i = 1:n_rings
    tmp=r(i)*exp(1j*((0:ntheta-1)*dtheta+mod(i,2)*dtheta/2));
    k   =   [k;real(tmp)' imag(tmp)'];
end
