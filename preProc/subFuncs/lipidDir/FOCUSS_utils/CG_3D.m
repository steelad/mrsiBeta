function [ x ] = CG_3D( A, W, y, lambda_, max_iter )

% try to modify Lustig's grad descent to minimize:
% || y - A*Wk*q ||^2 + lambda * ||q||^2  
% where y is the undersampled k-space data
%       A is the undersampled DFT operator
%       Wk is the matrix of Focuss weights at kth Focuss iteration
%       and q is the least squares soln to   A*Wk*q = y


x0 = zeros(size(W));


X = get_XFM(x0, 0, A, W);

g_o = grad_CG( X, x0, y, W, lambda_, A );
s = -g_o;
x = x0;
g_new = g_o;
dot_new = dot(g_new(:),g_new(:));
t0 = 1;

if nargin < 5
    max_iter = 20;
end

for k = 1:max_iter

    % line search
    alfa_ = t0;
    
    % precompute transforms
    [X,S] = get_XFM(x, s, A, W);
    
    func_LHS = cost_function_CG(X + alfa_*S, x + alfa_*s, y, lambda_);
    func_RHS1 = cost_function_CG(X, x, y, lambda_);
    func_RHS2 = .01 * abs( dot(g_new(:),s(:)) );
    
    LS_iter = 0;
    
    while func_LHS > func_RHS1 + alfa_ * func_RHS2  && LS_iter < 50
        alfa_ = .6 * alfa_;
        func_LHS = cost_function_CG(X + alfa_*S, x + alfa_*s, y, lambda_);
        LS_iter = LS_iter + 1;
    end
    
    
    if LS_iter > 2
        t0 = t0 * .6;
    elseif LS_iter < 1
		t0 = t0 / .6;
    end
    
    if mod(k,10) == 0;        disp(['Objective:  ', num2str(func_LHS)]);    end

    dot_o = dot_new;
   
    x = x + alfa_ * s;
    g_new = grad_CG( X + alfa_*S, x, y, W, lambda_, A );
    
    dot_new = dot(g_new(:),g_new(:));
    
    beta_ = dot_new / dot_o;
    s = -g_new + beta_ * s; 
   
end

disp(['Objective:  ', num2str(func_LHS)]);

end


% -------------------------------------------------------------------------
function [ g ] = grad_CG( X, x, y, W, lambda_, A )

g = 2 * ( lambda_*x + conj(W) .* ( A'*( X - y ) ) ); 

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function [ cost ] = cost_function_CG( X, x, y, lambda_ )

temp = X - y;
cost = norm(temp(:))^2;
cost = cost + lambda_ * norm(x(:))^2;

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function [ X, S ] = get_XFM( x, s, A, W )

    X = A*( x.*W );
    
    S = 0;
    if s ~= 0
       S = A*( s.*W );
    end
end
% -------------------------------------------------------------------------