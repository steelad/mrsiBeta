%***************************************************************************
%                                read_database.M
%***************************************************************************
% PURPOSE:  read database written in MAT format
%***************************************************************************
% CALL:     [met_prof] = read_database(dbase_file,met_names,table,ndp,)
%***************************************************************************
% INPUT:    dbase_file  =  name of the database file
%           met_names   =  metabolite names (ex: NAA,Lip1, etc)   
%           table       =  filename of the table in the directory tables
%           ndp         =  number of points 
%
% OUTPUT:   met_prof    =  metabolite profiles (each column corresponds to one 
%                          metabolite profile in the basis set) 
%***************************************************************************
% COMMENTS: 
%***************************************************************************

[met_prof] = read_database(dbase_file,ndp,met_names,table)

fid = fopen(table,'r');
refnames = {};
for i=1:noc
refnames(i) = cellstr(fscanf(fid,'%s',[1,1]));
end 
fclose(fid);

