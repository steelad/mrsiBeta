function [F,J] = ...
    funjaco(alf,dummy,stof,t,demptype,gausstype,freqtype,etype,...
            amplref,dempref,gaussref,freqref,phasref,eref,filter_info,...
            demp,gauss,freq,e,loline,galine,voline,t0,phasedistort,equalph,...
            data,filter_information,baseline_information,step,plotting)
% Computes current function value and jacobian at parameter 'alf',
% to be passed to the nonlinear least-squares solver

% dimensions
[nrp,nrstofs] = size(stof);
nrpf          = length(data);
nrvar         = length(alf);

j = sqrt(-1);

% read baseline variables, if baseline is used

if (~isempty(baseline_information))
  baseline_boolean = 1;
  B = baseline_information.B;
  D  = baseline_information.D;
  nosplpars = baseline_information.nosplpars;
  sord   =  baseline_information.sord;
  lambda = baseline_information.lambda;
else 
  baseline_boolean = 0;
end

[ampl,demp,gauss,freq,phas,e] = ...
    update(alf,nrstofs,demptype,gausstype,freqtype,etype,...
           amplref,dempref,gaussref,freqref,phasref,eref,demp,...
           gauss,freq,e,loline,galine,voline,phasedistort,equalph);


for i = 1:(nrstofs)
  A(:,i) = stof(:,i).*(exp((-demp(i)-gauss(i).*t+j*(freq(i)+e(i).*t)).*t)).';
end

if(~equalph)
  
  if (baseline_boolean == 1)
    Atot = [A B]; 
  else 
    Atot = A;
  end 

  % filter
  
  Atemp = filtering(Atot,filter_info,t);
  nrpf = size(Atemp,1);

  nrlinpars = nrstofs;

  if (baseline_boolean == 1)
    Atemp = [Atemp; zeros(size(D,1),nrstofs) lambda*D];
    nrlinpars = nrstofs + nosplpars;

  end

  [Q,R] = qr(Atemp); 
  Q1 = Q(:,1:nrlinpars); Q2 = Q(:,nrlinpars+1:end); 
  R  = R(1:nrlinpars,1:nrlinpars);
  zlin  = R\(Q1'*data);
  zampl = zlin(1:nrstofs);
  c     = zlin(nrstofs+1:end);
  
  % function value

  phi = Q2'*data;

  F(1:2:2*length(phi)-1) = real(phi);
  F(2:2:2*length(phi))   = imag(phi);
  F = F(:);

else
  
  zampl = ampl*exp(j*phas(1));
  recon = A*zampl(:);
  
  % filter  
  filrecon = filtering(recon,filter_info,t);
  nrpf = size(filrecon,1);

  if (baseline_boolean == 1)

    Atemp = [filtering(B,filter_info,t); lambda*D];

    [Q,R] = qr(Atemp); 
    Q1 = Q(:,1:nosplpars); Q2 = Q(:,nosplpars+1:end); 
    R  = R(1:nosplpars,1:nosplpars);
    filrecon = [filrecon;zeros(size(baseline_information.D,1),1)];
    c = R\(Q1'*(data-filrecon));
  
  else 
    Q2 = eye(nrpf);
  end


  % function value

  phi = Q2'*(data-filrecon);

  F(1:2:2*length(phi)-1) = real(phi);
  F(2:2:2*length(phi))   = imag(phi);
  F = F(:);
 
end


%
% computation of the Jacobian
%

dphi(1:nrp,1:nrvar) = zeros(nrp,nrvar);

% derivative wrt amplitude

if (equalph)
    for i = 1:nrstofs
      col = amplref(i);
      dphi(:,col) = A(:,i)*exp(j*phas(i));
    end
end

% derivative wrt damping

if (loline | voline)
    for i = 1:nrstofs
        select = demptype(i,1);
        
        if (select == 0 | select == 1)
            col = dempref(i);
            dphi(:,col) = dphi(:,col)-A(:,i)*zampl(i).*t';
        end
    end
end

%derivative wrt gaussian part damping

if (galine | voline)
    for i = 1:nrstofs
        select = gausstype(i,1);
        if (select == 0 | select == 1)
            col = gaussref(i);
            dphi(:,col) = dphi(:,col)-A(:,i)*zampl(i).*(t.^2)';
        end
    end
end


% derivative wrt frequency

for i = 1:nrstofs
    select = freqtype(i,1);
    
    if (select == 0 | select == 1)
        col = freqref(i);
        dphi(:,col) = dphi(:,col)+j*A(:,i)*zampl(i).*t';
    end
end


% derivative wrt phase

if (equalph)
    for i = 1:nrstofs
      col = phasref(i);
      dphi(:,col) = dphi(:,col)+j*A(:,i)*zampl(i);
    end
end


% derivative wrt e - residual eddy current effect

if (phasedistort)
    for i = 1:nrstofs
        select = etype(i,1);
        
        if (select == 0 | select == 1)
            col = eref(i);
            dphi(:,col) = dphi(:,col)+j*A(:,i)*zampl(i).*(t.^2)';
        end
    end
end

dphi2 = filtering(dphi,filter_info,t);

if  (baseline_boolean == 1)

  Jtemp = -Q2'*[dphi2; zeros(size(D,1),nrvar)];

else
  
  Jtemp = -Q2'*dphi2;

end
               
J(1:2:2*length(phi)-1,:) = real(Jtemp);
J(2:2:2*length(phi),:)   = imag(Jtemp);


% plot current fit
[xpts,xpts2]  = absi(nrpf,step); 
if (isfield(filter_info,'fir_h'))
  h = filter_info.fir_h;
else 
  h = 1;
end
[hr,w]     = freqz(h,1,nrpf,'whole');
hr(nrpf+1) = hr(1);

if (plotting == 1)
  figure(55);
  recon = A*zampl(:);
  recon = filtering(recon,filter_info,t);
  [cmp] = preplot(recon(1:nrpf).',hr,0,t0,xpts);
  [inp] = preplot(data(1:nrpf).',hr,0,t0,xpts); 
  plot(xpts2,inp,'b'), hold on
  plot(xpts2,cmp,'g')
  if (baseline_boolean == 1)
    baseline = filtering(B*c,filter_info,t);  
    [bsp] = preplot(baseline.',hr,0,t0,xpts);
    plot(xpts2,bsp,'black')
    plot(xpts2,bsp+cmp,'r')
  end
  set(gca, 'xdir', 'reverse')
  xlabel('a.u.')
  v = axis;
  % modified by JB
  %axis([0 4 v(3) v(4)]);  
  if baseline_boolean == 1
      legend('Original','Modeled metabolites','Baseline','Modeled signal (Green+black)')
  else
      legend('Original','Modeled metabolites')
  end
  grid;
  drawnow
  hold off
end



