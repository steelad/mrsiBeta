Matlab Toolbox for Dimensionality Reduction (v0.1b)
===================================================


Information
-------------------------
Author: Laurens van der Maaten
Affiliation: MICC-IKAT, Maastricht University, The Netherlands
Contact: l.vandermaaten@micc.unimaas.nl
Realease date: March 13th 2007
Version: 0.1b


Installation
------------
Copy the drtoolbox/ folder into the $MATLAB_DIR/toolbox directory (where $MATLAB_DIR indicates your Matlab installation directory). Start Matlab and select Set path... from the File menu. Click the Add with subfolders... button, select the folder $MATLAB_DIR/toolbox/drtoolbox in the file dialog, and press Open. Subsequently, press the Save button in order to save your changes to the Matlab search path. The toolbox is now installed. 
Some of the functions use MEX-files. Precompiled versions of these MEX-files are distributed with this release, but the compiled version for your platform might be missing. In order to compile all MEX-files, type CD $MATLAB_DIR/toolbox/drtoolbox in your Matlab prompt, and execute the function MEXALL. 


Features
-------------------------
This Matlab toolbox implements the following algorithms:

 - Principal Component Analysis ('PCA')
 - Linear Discriminant Analysis ('LDA')
 - Independent Component Analysis ('ICA')
 - Isomap ('Isomap')
 - Landmark Isomap ('LandmarkIsomap')
 - Locally Linear Embedding ('LLE')
 - Locally Linear Coordination ('LLC')
 - Laplacian Eigenmaps ('Laplacian')
 - Hessian LLE ('HessianLLE')
 - Local Tangent Space Alignment ('LTSA')
 - Diffusion maps ('DiffusionMaps')
 - Kernel PCA ('KernelPCA')
 - Generalized Discriminant Analysis ('KernelLDA')
 - Stochastic Neighbor Embedding ('SNE')
 - Neighborhood Preserving Embedding ('NPE')
 - Linearity Preserving Projection ('LPP')
 - Stochastic Proximity Embedding ('SPE')
 - Maximum Variance Unfolding ('MVU')
 - Autoencoders using RBM pretraining ('AutoEncoder')
 - Autoencoders using evolutionary optimization ('AutoEncoderEA')


Usage
-------------------------
Basically, you only need one function: mappedX = compute_mapping(X, technique, no_dims);
Try executing the following code:

	[X, labels] = generate_data('swiss');
	mappedX = compute_mapping(X, 'LLE', 2);
	scatter(mappedX(:,1), mappedX(:,2), 5, labels);

It will create a Swiss roll dataset, run LTSA on this dataset, and plot the result. Possible options for the techniques are 'PCA', 'LDA', 'ICA', 'LandmarkIsomap', 'Isomap', 'LLE', 'LLC', 'Laplacian', 'HLLE', 'LTSA', 'DiffusionMaps', 'KernelPCA', 'GDA', 'SNE', 'SPE', 'LPP', 'NPE', 'MVU', 'AutoEncoder', and 'AutoEncoderEA'. The COMPUTE_MAPPING function can also work directly on PRTools datasets (http://prtools.org). For more information, type HELP COMPUTE_MAPPING in your Matlab prompt.

Other functions that are useful are the GENERATE_DATA function and the TRANSFORM_SAMPLE_EST function. The generate_data function provides you with a number of artificial datasets to test the techniques. The TRANSFORM_SAMPLE_EST function allows you to perform an out-of-sample extension using a simple estimation technique.


Pitfalls
-------------------------
When you run certain code, you might receive an error that a certain file is missing. This is because in some parts of the code, MEX-functions are used. I provide a number of precompiled versions of these MEX-functions in the toolbox. However, the MEX-file for your platform might be missing. To fix this, type in your Matlab:

	mexall

Now you have compiled versions of the MEX-files as well. This fix also solves slow execution of the shortest path computations in Isomap.
If you encounter an error considering CSDP while running the FastMVU-algorithm (e.g., when running on a PowerPC Mac), the binary of CSDP for your platform is missing. If so, please obtain a binary distribution of CSDP from https://projects.coin-or.org/Csdp/ and place it in the drtoolbox/techniques directory. Make sure it has the right name for your platform (csdp.exe for Windows, csdpmac for Mac OS X, and csdplinux for Linux).

Many methods for dimensionality reduction perform spectral analyses of sparse matrices. You might think that eigenanalysis is a well-studied problem that can easily be solved. However, eigenanalysis of large matrices turns out to be tedious. The toolbox allows you to use two different methods for eigenanalysis:

	- The original Matlab functions (based on Arnoldi methods)
	- The JDQR functions (based on Jacobi-Davidson methods)

For problems up to 10,000 datapoints, we recommend using the 'Matlab' setting. For larger problems, switching to 'JDQR' is often worth trying.


Disclaimer
-------------------------
(C) Laurens van der Maaten, Maastricht University, 2007
You are free to use, modify, or redistribute this code in any way you want. If you do so, I would appreciate it if you mention the name of the author.

Parts of the code were taken from other authors, but often I made numerous improvements and modifications. A list of files in which I use source code from other authors is given below:
 - minimize.m: C.E. Rasmussen
 - hlle.m, mgs.m: C. Grimes and D. Donoho
 - autoencoder.m, backprop.m, rbm.m, rbmhidlinear.m, cg_update.m: G. Hinton and R. Salakhutdinov
 - dijk.m: M.G. Kay
 - dijkstra.cpp: J. Boyer
 - L2_distance.m: R. Bunschoten
 - jdqr.m, jdqz.m: G. Sleijpen
 - components.m: J. Gilbert
 - hillclimber2c.m, fastmvu.m, computegr.c, csdp.m, mexCCACollectData2.c, writesdpa.m, sparse_nn.m, readsol.m, sdecca2.m, hill_obj.m: K. Weinberger
 - llc.m, infermfa.m, mppca.m: Y. Teh


Contact
-------------------------
If you have any bugs, questions, suggestions, or modifications, please contact me:

	l.vandermaaten@micc.unimaas.nl

