function [dbase,abs_scl_factor,name_metabolite,index] = ...
    select_from_dbase(dbase0,abs_scl_factor0,name_metabolite0,...
                      dbase_selection)
% select some specified components from the database

nrp = size(dbase0,1);
nrstofs = 0;

for i=1:length(dbase_selection)

  index(i) = find(strcmp(dbase_selection{i},name_metabolite0) == 1);
  nrstofs = nrstofs + 1;
  name_metabolite{nrstofs} = name_metabolite0{index(i)};
  dbase(1:nrp,nrstofs)     = dbase0(1:nrp,index(i));
  abs_scl_factor(nrstofs,1)  = abs_scl_factor0(index(i));
  
end

% end select_from_dbase