function [ X ] = focuss_3D( X0, y, A, lambda_, cg_iter, iter_lim, p, Xtrue, lipid_mask)

% X0 : initial guess 
% y : undersampled k-space data
% A : DFT operator
% lambda_ : regularization parameter
% cg_iter : no of conj grad iterations to use in CG_3D
% iter_lim : no of outer Focuss iterations
% p : norm to impose on recon  (p=1 for L1, p<1 for non-convex recon) 
% Xtrue : fully sampled, true image (to print RMSE)


if nargin < 4
    lambda_ = 0;
end

if nargin < 5
    cg_iter = 50;
end

if nargin < 6
    iter_lim = 10;
end

if nargin < 7
    p = 1;
end

if nargin < 8
    Xtrue = 0;
end

if nargin < 9 
    lipid_mask = 0;
else
    lipid_mask = repmat(lipid_mask,[1,1,size(X0,3)]);
end


L = size(y,4);  % number of voxels to reconstruct jointly
X = X0;


W = sqrt( sum(X.*conj(X),4) ) .^ (1-p/2);


k = 1;

while 1
  
    disp(['-------- Iteration: ', num2str(k),' --------'])
     
    for m = 1:L
        q = CG_3D( A, W, y(:,:,:,m), lambda_, cg_iter );
        X(:,:,:,m) = (W.*q);
    end
       

    
    if Xtrue == 0
        figure(1), imagesc(sum(abs(X),3)), axis image, colormap jet,  drawnow
    else
        rmse_n = 100*norm( X(:).*lipid_mask(:)- Xtrue(:) ) / norm( Xtrue(:) );
        figure(1), imagesc(sum(abs(X),3)), axis image, colormap jet,  title(['Focuss RMSE: ', num2str(rmse_n), ' percent']), drawnow
    end
    
    
    W = sqrt( sum(X.*conj(X),4) ).^(1-p/2) ;
          

    if k >= iter_lim;        break;    end
      
    k = k+1; 
end


disp('-------- End of reconstruction --------')

if Xtrue ~= 0
    disp(['RMSE: ', num2str(100*norm(X(:)- Xtrue(:) ) / norm( Xtrue(:)) ), ' percent'])
end


end

