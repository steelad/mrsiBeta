function T1Spec = resample_spectroscopy_slice( T1_file, T1_BetMask_file, Spec_file, Save_folder )
%
% T1Spec = resample_spectroscopy_slice( T1_file, T1_BetMask_file, Spec_file )
% T1Spec = resample_spectroscopy_slice( T1_file, T1_BetMask_file, Spec_file, Save_folder )
%
% Takes a T1 image, the corresponding binary BET mask, and a spectroscopy slice.
% Then for each metabolite, resample the spectroscopy slice in T1 space.
% Output is a cell-array, with one T1 volume for each metabolite.
%
% If Save_folder is specified, then each metabolite volume is saved in that folder.
%
% JH

    % Load both images
    disp('Loading images...');
    T1  = load_untouch_nii(T1_file);
    T1m = load_untouch_nii(T1_BetMask_file);
    Sp  = load_untouch_nii(Spec_file);
    
    % coordinates of brain voxels
    disp('Working out resampling...');
    [T1ind,~,T1xyz] = myfind(T1m);
    
    % corresponding indices in spectroscopy slice
    Spijk = nifti_xyz2ijk( Sp, T1xyz );
    Spsiz = size(Sp.img);
    try
        Nmeta = Spsiz(4); % number of metabolites
    catch
        warning('++ One metabolite in the map. Coercing map-values size')
        Spsiz([3 4]) = [1 1]
        Nmeta = 1;
    end
    Spnum = prod(Spsiz(1:3));
    
    % filter indices that are within the slice
    Spind = mysub2ind( Spsiz, Spijk );
    valid = (Spind >= 1) & (Spind <= Spnum);
    Spind = Spind(valid);
    T1ind = T1ind(valid);
    
    % create output volumes
    T1_zero = T1;
    T1_zero.img = zeros( size(T1.img), class(T1.img) );
    T1Spec = cell(1,Nmeta);
    
    for i = 1:Nmeta
        fprintf('\tMetabolite #%02d...\n',i);
        Ti = T1_zero;
        Spi = Sp.img(:,:,:,i);
        Ti.img = single(Ti.img);
        Ti.hdr.dime.datatype = 64;
        Ti.img(T1ind) = Spi(Spind);
        T1Spec{i} = Ti;
        if nargin > 3
            if i == 1
                mkdir(Save_folder)
            end
            save_untouch_nii( Ti,fullfile(Save_folder,sprintf('T1Spec_metabolite%02d.nii',i)));
        end
    end
    
end

function [ind,ijk,xyz] = myfind( nii )
    assert( ndims(nii.img)==3 || size(nii.img,4)==1, 'Input image should be 3D.' );
    ind = find(nii.img);
    [i,j,k] = ind2sub( size(nii.img), ind );
    ijk = [i,j,k];
    xyz = nifti_ijk2xyz( nii, ijk );
end

function ind = mysub2ind( imsize, ijk )
    imsize = [1,cumprod(imsize)];
    ind = 1 + sum(bsxfun( @times, ijk-1, imsize(1:3) ),2);
end
