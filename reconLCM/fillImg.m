function new = fillImg(LUT,metab_map2,trans)
new  = zeros([size(LUT) size(metab_map2,3)]);
if trans == 1
        metab_map2 = flipud(metab_map2);
    elseif trans == 2
        metab_map2 = rot90(metab_map2,-90);
        metab_map2 = circshift(metab_map2,-1,2);
%         metab_map2 = circshift(metab_map2,-1,1);
end
for j = 1:size(metab_map2,3)
    disp(['++ Starting slice ' num2str(j)])
<<<<<<< HEAD
    if trans == 1
        %metab_map2 = permute(metab_map2,[2 1 3]);
        metab_map2 = flipud(metab_map2);
    end
=======
    
>>>>>>> dev
    B = metab_map2(:,:,j);
    B = B(:);
    ind = 1:length(B);
     ind = ind(B~=0); B = B(B~=0);
    %%
    temp = zeros(size(LUT));
    for i = 1:length(B);
        temp(LUT == ind(i)) = B(i);
    end
    %%
    new(:,:,:,j) = temp;
    clear temp ind
end
disp(['++ output size: ' num2str(size(new))])
%%
