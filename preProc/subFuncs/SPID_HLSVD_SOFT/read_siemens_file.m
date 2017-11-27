function [signal, ndp, step, nfids]= read_siemens_file(file2read)


step  = 1;
nfids = 1;
fin   = fopen(file2read,'r','ieee-le');
ndp   = 1024;
[A,count] = fread(fin,[2 ndp],'float');
signal = A(1,1:ndp) + i * A(2,1:ndp);
fclose(fin)
