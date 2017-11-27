%***************************************************************************
%                                HLSVDPRO.M
%***************************************************************************
% PURPOSE:  Estimation algorithm with LS solution, in U space, for the  
%           nonlinear parameters and LS solution for the linear
%           parameters.
%           Computation of the singular values by the Lanczos
%           algorithm with partial reorthogonalization (PRO).
%***************************************************************************
% CALL:     [freq,damp,ampl,phas]=HLSVDPRO(vhsvd1,vhsvd2,Kest,fsampl,t,M)
%***************************************************************************
% INPUT:    vhsvd1 -- chosen data vector for the 1st system of HSVD algorithm
%           vhsvd2 -- chosen data vector for the 2nd system of HSVD algorithm
%           Kest   -- estimated number of damped exponentials (model order)
%           fsampl -- sample frequency
%           t      -- sampled time sequence for vhsvd
%           M      -- number of matrix columns
% OUTPUT:   freq   -- estimated frequency      row vector
%           damp   -- estimated damping factor row vector
%           ampl   -- estimated amplitude      row vector
%           phas   -- estimated phase          row vector
%***************************************************************************

function [freq,damp,ampl,phas,SV]=HLSVDPRO(vhsvd1,vhsvd2,Kest,fsampl,t,M)
%global harrayflops1

%global harrayflops2
%flops(0)

%-------------------------------initialization------------------------------
N=length(vhsvd1); L=N+1-M;
vhsvd1=vhsvd1(:);

 colhmat=vhsvd1(N-L+1:N);
  rowhmat=vhsvd1(N-L+1:-1:1);
%flops(0);
[u,s,v,info1,info2]=lansvdMNT(colhmat,rowhmat,Kest,'L');
%fl =flops;

%-------------------------Truncation to the given order---------------------
u1=u(:,1:Kest);
s1=s(1:Kest,1:Kest);
v1=v(:,1:Kest);

%f1=flops;
%flops(0)

%------------------------LS solution for the 1st system---------------------
up=u1;
upt=up(1:L-1,:);           % top part of matrix [up] without bottum line 
upb=up(2:L,:);             % bottom part of matrix [up] without top line
q=pinv(upt)*upb;           % LS solution 

z=eig(q).';
[s, I]=sort(imag(log(z))); % arranges the frequencies in ascending order
z=z(I);
a=fsampl/(2*pi);
for i=1:Kest
    freq(i)=imag(log(z(i))*a);
    damp(i)=-real(log(z(i))*fsampl);
end
%---------------------LS solution for the 2nd system------------------------

[ampl,phas]=linpara(vhsvd2,Kest,freq,damp,t);

%f2=flops;
%harrayflops1=[harrayflops1;f1];
%harrayflops2=[harrayflops2;f2];

if nargout < 5
    return
else
    SV = s;
end
