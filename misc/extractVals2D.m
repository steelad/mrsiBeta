function vals = extractVals2D(metab_map2,mask,metab)
%%
temp = metab_map2(:,:,metab);
mV = unique(mask); mV = mV(mV ~= 0);
vals = cell(1,length(mV));

for i = 1:length(mV)
    vals{i} = temp(mask == mV(i)); vals{i} = vals{i}(vals{i}~=0);
end
return