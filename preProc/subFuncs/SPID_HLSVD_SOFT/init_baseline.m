function baseline_information = init_baseline(xlow,xhi,ldata,ndx,bdeg,...
                                              sord,li,lambda)

nosplpars = ndx+bdeg;             % number of spline pars

xx = [1:ldata]';                  % equidistant abscissas
xl = xlow ; xr = xhi;  
B  = bsplinep(xx,xl,xr,ndx,bdeg); % B-spline basis
B(:,1:bdeg) = zeros(ldata,bdeg); 
B(:,end-bdeg+1:end) = zeros(ldata,bdeg);

% store (in baseline_information) all values to be passed to the solver

baseline_information.B = ifft(B);
baseline_information.D = [diff(eye(nosplpars),sord); ...
                          -1 zeros(1,nosplpars-2) 1];
baseline_information.nosplpars = nosplpars;
baseline_information.sord = sord;
baseline_information.lambda = lambda;

if lambda < 0,

  [U,V,X,alpha,beta] = gsvd(B,baseline_information.D,0);
  baseline_information.alpha = diag(alpha'*alpha);
  baseline_information.beta  = diag( beta'*beta );
  
end
