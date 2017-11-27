function T1_img = T1_import_wrapper(T1_file)
%%
% wrapper for nifit toolbox import t1, keeps it in radiological convention
% and rotates for mrsi
%%

T1_file = load_nii(T1_file); 
T1_img = T1_file.img;
T1_img = rot90(T1_img);
T1_img = fliplr(T1_img);
