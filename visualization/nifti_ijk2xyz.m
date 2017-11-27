function [xyz,T] = nifti_ijk2xyz( nii_struct, ijk )

    assert( size(ijk,2) == 3, 'ijk should be nx3.' );

    % get orientation matrix
    T = nifti_orientation( nii_struct );
    
    % apply transformation
    N   = size(ijk,1);
    ijk = horzcat( ijk-1, ones(N,1) ); % Matlab indices begin at 1
    xyz = ijk * T(1:3,:)';
    
end
