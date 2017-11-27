function [xyz,T] = nifti_mask2coord( nii_struct, use_untouch )

    if nargin < 2, use_untouch = true; end

    % if input is a filename, load it
    if ischar( nii_struct )
        if use_untouch
            nii_struct = load_untouch_nii(nii_struct);
        else
            nii_struct = load_nii(nii_struct);
        end
    end
    
    % make sure that it's a 3D (not 4D) volume
    switch ndims(nii_struct.img)
        
        case 3
            % get ijk indices of non-zero voxels
            [i,j,k] = ind2sub( size(nii_struct.img), find(nii_struct.img) );
            
        case 4
            warning('Selecting non-zero voxels across 4th dimension to find 3D coordinates.');
            
            img_size = size(nii_struct.img);
            img_size = img_size(1:3);
            
            [i,j,k] = ind2sub( img_size, find(max(nii_struct.img,[],4)) );
            
        otherwise
            error('Unsupported data dimensions.');
    end
    
    % apply transformation
    [xyz,T] = nifti_ijk2xyz( nii_struct, [i(:),j(:),k(:)] );
    
end
