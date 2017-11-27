function T = nifti_orientation( nii_struct, method )

    % extract relevant data
    dat    = nii_struct.hdr.hist;
    qfac   = nii_struct.hdr.dime.pixdim(1);
    pixdim = nii_struct.hdr.dime.pixdim([2 3 4]);
    
    % method of conversion
    if nargin < 2
        if dat.qform_code == 0 && dat.sform_code == 0
            method = 1; % homothetic
        elseif dat.qform_code > 0
            method = 2; % orthogonal
        else
            method = 3; % affine
        end
    end
    
    % compute transformation matrix
    switch method
        
        case 1
            
            T = horzcat(  diag(pixdim), [0;0;0]  );
            
        case 2
            
            % quaternion
            b = dat.quatern_b;
            c = dat.quatern_c;
            d = dat.quatern_d;
            a = 1 - b*b - c*c - d*d;
            
            % deal with precision errors
            if a < 0
                warning( 'Quaternion does not have unit norm (a^2=%g)', a );
                a = 0;
            else
                a = sqrt(a);
            end
            
            % offset
            ox = dat.qoffset_x;
            oy = dat.qoffset_y;
            oz = dat.qoffset_z;
            
            % build transformation matrix
            S = pixdim .* [1 1 qfac];
            R = [ 1 - 2*c*c - 2*d*d , 2*b*c - 2*d*a     , 2*b*d + 2*c*a ; ...
                  2*b*c + 2*d*a     , 1 - 2*b*b - 2*d*d , 2*c*d - 2*b*a ; ...
                  2*b*d - 2*c*a     , 2*c*d + 2*b*a     , 1 - 2*c*c - 2*b*b ];
            T = horzcat(  bsxfun(@times,R,S), [ox;oy;oz]  );
            
        case 3
            
            T = [ dat.srow_x; dat.srow_y; dat.srow_z ];
        
    end

    % append the last row to make T a square matrix
    T = [T; 0 0 0 1];
    
end
