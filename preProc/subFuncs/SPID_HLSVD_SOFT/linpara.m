%***************************************************************************
%                                 linpara.m
%***************************************************************************
% PURPOSE:  Estimate linear paramenters of sum of exponentially damped 
%           sinusoids using LS technique.
%***************************************************************************
% CALL:     [ampl,phas]=linpara(vlin,Klin,freqlin,damplin,tlin)
%***************************************************************************
% INPUt:    vlin   -- chosen data vector for estimating linear parameters
%           Klin   -- estimated number of damped exponentials (model order)
%           freqlin-- estimated frequency row vector
%           damplin-- estimated damping factor row vector
%           tlin   -- sampled time sequence for vlin
% OUTPUT:   ampl   -- estimated amplitude      row vector
%           phas   -- estimated phase          row vector
%***************************************************************************

function [ampl,phas]=linpara(vlin,Klin,freqlin,damplin,tlin)

%-------------------LS solution for the 2nd system--------------------------
b=tlin'*(-damplin+sqrt(-1)*2*pi.*freqlin);
x=exp(b);

% weighting
w=[];
for i=1:Klin
  w=[w norm(x(:,i))];
end
w=diag(w);
if cond(w) > 10^14
  disp('my warning in linpara(1).m: the cond(w) > 10^14. ')
end
w=inv(w);
x=x*w;

[u,s,v]=svd(x,0);
u1=u(:,1:Klin);
s1=s(1:Klin,1:Klin);
v1=v(:,1:Klin);

if cond(s1) > 10^14
  disp('my warning in linpara.m(2): the cond(s1) > 10^14. ')
end



c=w*v1*inv(s1)*u1'*vlin(:);           % with weighting
%c=v1*pinv(s1)*u1'*vhsvd2(:);         % if without weighting

ampl=abs(c).';
y=imag(c); x=real(c);
%----------------------------adjustment of phase----------------------------
for i=1:Klin,
  if x(i)==0,
    if y(i)>=0,
      phas(i)=90;
    else
      phas(i)=270;
    end
  else,
    phas(i)=atan(y(i)/x(i))*180/pi;
  end
end

% the following makes phas correct and positive : from 0 to 360- 
for i=1:Klin,
  if y(i)<0 & x(i)>0;
    phas(i)=phas(i)+360;
  end
  if x(i)<0;                      % necessary when x<0 and y>0.
    phas(i)=phas(i)+180;
  end
end
% to change phas from 0 -- 360 to -180 -- 180-, use phaseadj or 
% delete the above lines and change them to if x<0 & y>0, phas+180, end
