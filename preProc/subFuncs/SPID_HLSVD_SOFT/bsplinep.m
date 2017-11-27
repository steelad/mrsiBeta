function B= bsplinep(x,xl,xr,ndx,bdeg)
% x: input value for which to compute basis function  
% xl: leftmost point
% xr: rightmost point
% ndx: number of intervals
% bdeg: degree of basis functions of B-spline
  
dx = (xr-xl)/ndx;% knot interval length
t  = xl+dx*[-bdeg:ndx-1];% knot sequence
T  = (0*x+1)*t;% knot sequence repeated as many times as there are
               % input values
X  = x*(0*t+1);% input vector repeated as many times as there are knots
P  = (X-T)/dx;
B  = (T<=X)&(X<(T+dx)); 
r  = [2:length(t) 1];
for k = 1:bdeg
  B = (P.*B+(k+1-P).*B(:,r))/k;
end
