function T1sp = resampling()

    % Load both images
    disp('Loading images...');
    FOLDER='/Users/jhadida/Data/temp/adam';
    IMAGES= { 'TUM_474_highres.nii.gz', '_TUM_474_MAP2.nii' };
    
    T1 = load_untouch_nii(fullfile( FOLDER, IMAGES{1} ));
    Sp = load_untouch_nii(fullfile( FOLDER, IMAGES{2} ));
    
    % Binary mask, showing all voxels with non-zero concentrations across metabolites
    Nmetabolites = size(Sp.img,4);
    Spmask = max( Sp.img, [], 4 ) > 0;
    
    % Upsample the mask
    disp('Upsampling...');
    Spfac  = 20;
    Spmask = imresize( Spmask, Spfac, 'nearest' );
    [subi,subj] = ind2sub( size(Spmask), find(Spmask) ); % indices in the upsampled image
    subi = 1 + (subi-1) / Spfac; 
    subj = 1 + (subj-1) / Spfac; 
    
    % Find their coordinate in T1 space
    xyz = cell( 1, Spfac );
    kval = 1:(1/Spfac):2;
    kval = kval(1:end-1); % < 2
    for i = 1:Spfac
       xyz{i} = nifti_ijk2xyz( Sp, [subi,subj,kval(i)*ones(size(subi))] ); 
    end
    %xyz = nifti_ijk2xyz( Sp, [i,j,ones(size(i))] ); % set k=1 for all pixels
    
    % For each of these points, find the nearest voxel in the T1 image
    disp('Finding voxel indices in original images...');
    T1sub = nifti_xyz2ijk( T1, vertcat(xyz{:}) );
    T1ind = mysub2ind( size(T1.img), T1sub );
    Spsub = nifti_xyz2ijk( Sp, xyz{1} );
    Spind = mysub2ind( size(Sp.img), Spsub );
    
    assert( all(Spsub(:,3) == 1), 'Something went wrong.' ); % k should be 1
    
    % Group points in xyz based on the voxel the fall into in the T1 image
    disp('Grouping points by voxel...');
    T1vox = unique(T1ind);
    T1vox = T1vox((T1vox >= 1) & (T1vox <= numel(T1.img))); % only valid voxels indices
    T1nvox = numel(T1vox);
    T1group = cell(1,T1nvox);
    
    for i = 1:T1nvox
        T1group{i} = find(T1ind == T1vox(i)); % find the index of the xyz points that fall in voxel T1vox(i)
    end
    
    % Create an empty image like the T1 for each metabolite
    T1new = T1;
    T1new.img = zeros( size(T1.img), class(T1.img) );
    T1sp = cell(1,Nmetabolites);
    
    for i = 1:Nmetabolites
        
        fprintf('\t Metabolite #%d...\n',i);
        
        % empty image for metabolite i
        Ti = T1new;
        
        % spectroscopy image for metabolite i
        Spi = squeeze(Sp.img(:,:,1,i));
        Spv = Spi(Spind); % values corresponding to each point in xyz
        Spv = repmat(Spv,Spfac,1);
        
        % for each voxel in T1, assign the most common value amongst all the xyz point that fall in it
        for j = 1:T1nvox
            Ti.img(T1vox(j)) = mode(Spv(T1group{j}));
        end
        
        % save that image
        T1sp{i} = Ti;
        file = sprintf( 'T1_metabolite%02d.nii.gz', i );
        save_untouch_nii( Ti, fullfile(FOLDER,file) );
        
    end

end

function ind = mysub2ind( s, ijk )

    s = [1,cumprod(s)];
    ind = 1 + sum(bsxfun( @times, ijk-1, s(1:3) ),2);
    
    % filter the points that are outside
    %ind = ind( (ind >= 1) & (ind <= s(4)) );

end
