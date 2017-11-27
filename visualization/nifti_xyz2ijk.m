function [ijk,T] = nifti_xyz2ijk( nii_struct, xyz )

    assert( size(xyz,2) == 3, 'xyz should be nx3.' );

    % get orientation matrix
    T = nifti_orientation( nii_struct );
    
    % reverse transformation (this could be optimized but meh)
    N   = size(xyz,1);
    xyz = horzcat( xyz, ones(N,1) ); 
    ijk = round( 1 + xyz / (T') ); % Matlab indices begin at 1
    ijk = ijk(:,1:3); % don't call unique(.,'rows') to preserve mapping with xyz

end
