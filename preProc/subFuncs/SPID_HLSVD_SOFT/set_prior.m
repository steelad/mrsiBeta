function [prior_equal,equal_to] = set_prior(prior_equal_name,equal_to_name,...
                                                          name_metabolite)
% function that returns indeces of components with equal dampings
% in the database, given their names.

ii = 0; prior_equal = []; equal_to = [];
for i = 1:length(prior_equal_name)

  temp1 = find(strcmp(name_metabolite,prior_equal_name(i))==1);
  temp2 = find(strcmp(name_metabolite,equal_to_name(i))==1);
  
  if (~isempty(temp1) & ~isempty(temp2))
    
    ii = ii + 1;
    prior_equal(ii) = temp1;
    equal_to(ii)    = temp2;
    
  end

end