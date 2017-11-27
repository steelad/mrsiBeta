function [term] = diff_sum(current_data,data,neighbor_vector,a,max_values,min_values)

    term = 0;
    for j = 1:size(neighbor_vector,2)
         term = term + (abs(current_data(a)-data(neighbor_vector(j),a))/(max_values(a)-min_values(a)));
    end;

end