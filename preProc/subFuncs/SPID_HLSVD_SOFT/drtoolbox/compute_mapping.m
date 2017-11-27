function [mappedA, mapping] = compute_mapping(A, type, no_dims, varargin)
%COMPUTE_MAPPING Performs dimensionality reduction on a dataset
%
%   mappedA = compute_mapping(A, type)
%   mappedA = compute_mapping(A, type, no_dims)
%
% Performs a technique for dimensionality reduction on the data specified 
% in A, reducing data with a lower dimensionality in mappedA.
% The data on which dimensionality reduction is performed is given in A
% (rows correspond to observations, columns to dimensions). A may also be a
% (labeled or unlabeled) PRTools dataset.
% The type of dimensionality reduction used is specified by type. Possible
% values are 'PCA', 'SPCA', 'LDA', 'ICA', 'MDS', 'Isomap', 'LandmarkIsomap',
% 'LLE', 'Laplacian', 'HessianLLE', 'LTSA', 'DiffusionMaps', 'KernelPCA', 
% 'GDA', 'SNE', 'SPE', 'AutoEncoder', and 'AutoEncoderEA'. 
% The function returns the low-dimensional representation of the data in the 
% matrix mappedA. If A was a PRTools dataset, then mappedA is a PRTools 
% dataset as well. For some techniques, information on the mapping is 
% returned in the struct mapping.
% The variable no_dims specifies the number of dimensions in the embedded
% space (default = 2). For 'LDA' and 'GDA', the labels of the instances 
% should be specified in the first column of A. 
%
%   mappedA = compute_mapping(A, type, no_dims, parameters)
%   mappedA = compute_mapping(A, type, no_dims, parameters, eig_impl)
%
% Free parameters of the techniques can be defined as well (on the place of
% the dots). These parameters differ per technique, and are listed below.
% For techniques that perform spectral analysis of a sparse matrix, one can 
% also specify in eig_impl the eigenanalysis implementation that is used. 
% Possible values are 'Matlab' and 'JDQR' (default = 'Matlab'). We advice
% to use the 'Matlab' for datasets of with 10,000 or less datapoints; 
% for larger problems the 'JDQR' might prove to be more fruitful. 
% The free parameters for the techniques are listed below (the parameters 
% should be provided in this order):
%
%   PCA:            - none
%   LDA:            - none
%   ICA:            - none
%   MDS:            - none
%   Isomap:         - <int> k -> default = 12
%   LandmarkIsomap: - <int> k -> default = 12
%                   - <double> percentage -> default = 0.2
%   LLE:            - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   LLC:            - <int> k -> default = 12
%                   - <int> no_analyzers -> default = 20
%                   - <int> max_iterations -> default = 200
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   Laplacian:      - <int> k -> default = 12
%                   - <double> sigma -> default = 1.0
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   HessianLLE:     - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   LTSA:           - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   MVU:            - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   DiffusionMaps:  - <double> sigma -> default = 1.0
%   KernelPCA:      - <char[]> kernel -> {'linear', ['gauss'], 'poly', 'subsets'} 
%                   - kernel parameters: type HELP GRAM for info
%   GDA:            - <char[]> kernel -> {'linear', ['gauss'], 'poly', 'subsets'} 
%                   - kernel parameters: type HELP GRAM for info
%   SNE:            - <double> sigma -> default = 1.0
%   LPP:            - none
%   NPE:            - none
%   SPE:            - <char[]> type -> {['Global'], 'Local'}
%                   - if 'Local': <int> k -> default = 12
%   AutoEncoder:    - none
%   AutoEncoderEA:  - none
%
% In the parameter list above, {.., ..} indicates a list of options, and []
% indicates the default setting. The variable k indicates the number of 
% nearest neighbors in a neighborhood graph. The variable sigma indicates
% the variance of a Gaussian kernel.
% For more information on the techniques, we refer to the paper
% "Dimensionality Reduction: A Comparative Review" by L.J.P. van der
% Maaten, E.O. Postma, and H.J. van den Herik. The paper is available from
% http://www.cs.unimaas.nl/l.vandermaaten.

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.1b.
% The toolbox can be obtained from http://www.cs.unimaas.nl/l.vandermaaten
% You are free to use, change, or redistribute this code in any way you
% want. However, it is appreciated if you maintain the name of the original
% author.
%
% (C) Laurens van der Maaten
% Maastricht University, 2007

    % Check inputs
    if nargin < 2
        error('Function requires at least two inputs.');
    end
    if ~exist('no_dims', 'var')
        no_dims = 2;
    end
    if length(varargin) > 0 && strcmp(varargin{length(varargin)}, 'JDQR')
        eig_impl = 'JDQR';
        varargin(length(varargin)) = [];
    elseif length(varargin) > 0 && strcmp(varargin{length(varargin)}, 'Matlab')
        eig_impl = 'Matlab';
        varargin(length(varargin)) = [];
    else
        eig_impl = 'Matlab';
    end
    mapping = struct;
    
    % Handle PRTools dataset
    if strcmp(class(A), 'dataset')
        prtools = 1;
        AA = A;
        if ~strcmp(type, {'LDA', 'FDA', 'GDA', 'KernelLDA', 'KernelFDA'})
            A = A.data;
        else
            A = [A.labels A.data];
        end
    else 
        prtools = 0;
    end
    
    % Switch case
    switch type
        case 'Isomap'         
            % Compute Isomap mapping
			if isempty(varargin), mappedA = isomap(A, no_dims, 12);
            else mappedA = isomap(A, no_dims, varargin{1}); end
			
		case 'LandmarkIsomap'
			% Compute Landmark Isomap mapping
            if isempty(varargin), mappedA = landmark_isomap(A, no_dims, 12, 0.2);
			elseif length(varargin) == 1, mappedA = landmark_isomap(A, no_dims, varargin{1}, 0.2);
            elseif length(varargin) >  1, mappedA = landmark_isomap(A, no_dims, varargin{1}, varargin{2}); end
            
        case {'Laplacian', 'LaplacianEigenmaps'}
            % Compute Laplacian Eigenmaps-based mapping
            if isempty(varargin), mappedA = laplacian_eigen(A, no_dims, 12, 1, eig_impl);
			elseif length(varargin) == 1, mappedA = laplacian_eigen(A, no_dims, varargin{1}, 1, eig_impl);
            elseif length(varargin) > 1,  mappedA = laplacian_eigen(A, no_dims, varargin{1}, varargin{2}, eig_impl); end
            
        case {'HLLE', 'HessianLLE'}
            % Compute Hessian LLE mapping
			if isempty(varargin), mappedA = hlle(A, no_dims, 12, eig_impl);
            else mappedA = hlle(A, no_dims, varargin{1}, eig_impl); end
            
        case 'LLE'
            % Compute LLE mapping
			if isempty(varargin), mappedA = lle(A, no_dims, 12, eig_impl);
            else mappedA = lle(A, no_dims, varargin{1}, eig_impl); end
            
        case 'LLC'
            % Compute LLC mapping
			if isempty(varargin), mappedA = run_llc(A', no_dims, 12, 20, 200, eig_impl);
            elseif length(varargin) == 1, mappedA = run_llc(A', no_dims, varargin{1}, 20, 200, eig_impl);
            elseif length(varargin) == 2, mappedA = run_llc(A', no_dims, varargin{1}, varargin{2}, 200, eig_impl);
            else mappedA = run_llc(A', no_dims, varargin{1}, varargin{2}, varargin{3}, eig_impl); end
            mappedA = mappedA';
           
        case 'LTSA'
            % Compute LTSA mapping 
            if isempty(varargin), mappedA = ltsa(A, no_dims, 12, eig_impl);
            else mappedA = ltsa(A, no_dims, varargin{1}, eig_impl); end
            
        case {'MVU', 'FastMVU'}
            % Compute MVU mapping
            if isempty(varargin), mappedA = fastmvu(A, no_dims, 12, eig_impl);
            else mappedA = fastmvu(A, no_dims, varargin{1}, eig_impl); end
            
        case {'DM', 'DiffusionMaps'}
            % Compute diffusion maps mapping
			if isempty(varargin), mappedA = diffusion_maps(A, no_dims, 1);
            else mappedA = diffusion_maps(A, no_dims, varargin{1}); end
            
        case 'SPE'
            % Compute SPE mapping
            if isempty(varargin), mappedA = spe(A, no_dims, 'Global');
            elseif length(varargin) == 1, mappedA = spe(A, no_dims, varargin{1}); 
            elseif length(varargin) == 2, mappedA = spe(A, no_dims, varargin{1}, varargin{2});
            end
            
        case 'LPP'
            % Compute LPP mapping
            mappedA = lpp(A, no_dims);
            
        case 'NPE'
            % Compute NPE mapping
            mappedA = npe(A, no_dims);
            
        case 'SNE'
            % Compute SNE mapping
			if isempty(varargin), mappedA = sne(A, no_dims);
            else mappedA = sne(A, no_dims, varargin{1}); end

        case 'AutoEncoder'
            % Train autoencoder to map data
            [mappedA, mapping] = autoencoder(A, no_dims);
            
        case 'AutoEncoderEA'
            % Train autoencoder to map data
            [mappedA, mapping] = autoencoder_ea(A, no_dims);
            
        case 'ICA'
            % Compute ICA mapping
            mappedA = ica(A);
            
        case 'KernelPCA'
            % Apply kernel PCA with polynomial kernel
			if isempty(varargin), [mappedA, mapping] = kernel_pca(A, no_dims);
			else [mappedA, mapping] = kernel_pca(A, no_dims, varargin{:}); end
			
		case {'KernelLDA', 'KernelFDA', 'GDA'}
			% Apply GDA with Gaussian kernel
            if isempty(varargin), mappedA = gda(A(:,2:end), uint8(A(:,1)), no_dims);
            else mappedA = gda(A(:,2:end), uint8(A(:,1)), no_dims, varargin{:}); end
            
        case {'LDA', 'FDA'}
            % Construct labeled dataset
            [mappedA, mapping] = lda(A(:,2:end), uint16(A(:,1)), no_dims);
            
        case 'MDS'
            % Compute distance matrix
            mappedA = mds(A, no_dims);   
            
        case 'PCA'
            % Compute Karhunen-Loeve mapping
			[mappedA, mapping] = pca(A, no_dims);   
            
        otherwise
            error('Unknown dimensionality reduction technique.');
    end
    
    % JDQR makes empty figure; close it
    if strcmp(eig_impl, 'JDQR')
        close(gcf);
    end
    
    % Handle PRTools dataset
    if prtools == 1
        AA.data = mappedA;
        mappedA = AA;
    end
    