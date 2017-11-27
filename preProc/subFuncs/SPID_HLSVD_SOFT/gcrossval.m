function gcv = gcrossval(F,alpha,beta,lambda,nrp,nosplpars)
  
% Function that computes the Generalized Cross Validation criterion
% for the regularized minimization of AQSES, at a given value of 
% the regularization parameter lambda
%
% Numerator of GCV function equals the norm of the residual error.
%
% Denominator of GCV function is the square of the trace
% of (I - H(lambda)), where H(lambda) is the smoother matrix
% H(lambda) = B*(B'*B + lambda^2 * D'*D)^(-1)*B' 
%            (B = coef. matrix, D = penalty matrix)
% Used simplification:                            lambda^2*beta(i)^2
%                    denominator = sum_i  -------------------------------
%                                          alpha(i)^2 + lambda^2*beta(i)^2


% Diana Sima, KUL-ESAT, May 2005


% Compute numerator

misfit = norm(F);

% Compute denominator

trace = ...%nrp - nosplpars + ...
        sum(lambda^2 * beta ./ (alpha + lambda^2 * beta));

% Compute GCV value

gcv = misfit / trace;

