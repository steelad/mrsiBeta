function [hits,missess] = find_neighbors(group1,group2,current,neighbors)

    total_data = [group1;group2;];
    current_data = total_data(current,:);
    
    max_values = max(total_data);
    min_values = min(total_data);
    
    row_number1 = size(group1,1);
    row_number2 = size(group2,1);
    
    group1_minus = [];
    for i = 1:size(group1,1)
        row = [];
        for j = 1:size(group1,2)
            row = [row abs(group1(i,j) - current_data(j))];
        end;
        group1_minus = [group1_minus; row;];        
    end
    
    group2_minus = [];
    for i = 1:size(group2,1)
        row = [];
        for j = 1:size(group2,2)
            row = [row abs(group2(i,j) - current_data(j))];
        end;
        group2_minus = [group2_minus; row; ];       
    end
    
    group1_minus_norm = [];
    for i = 1:size(group1_minus,2)
        col = group1_minus(:,i) / (max_values(i)-min_values(i));
        group1_minus_norm = [group1_minus_norm col];
    end

    group2_minus_norm = [];
    for i = 1:size(group2_minus,2)
        col = group2_minus(:,i) / (max_values(i)-min_values(i));
        group2_minus_norm = [group2_minus_norm col];
    end
    
    group1_diff = sum(group1_minus_norm');
    group2_diff = sum(group2_minus_norm');
        
    if (current > size(group1,1))
        [missess_total_sorted,missess_total] = sort(group1_diff);
        [hits_total_sorted,hits_total] = sort(group2_diff);
        missess = missess_total(1:neighbors);
        current_group2 = current - size(group1,1);
        hits = [];
                
        for i = 1:neighbors + 1
            if current_group2 == hits_total(i)
            else
                hits = [hits hits_total(i)];
            end;
        end;
        
    else
        [missess_total_sorted,missess_total] = sort(group2_diff);
        [hits_total_sorted,hits_total] = sort(group1_diff);
        missess = missess_total(1:neighbors);
        hits = [];
                
        for i = 1:neighbors + 1
            if current == hits_total(i)
            else
                hits = [hits hits_total(i)];
            end;
        end;
        
    end
    
    if size(hits) > neighbors 
        hits = hits(1:neighbors)
    
    end
    
end