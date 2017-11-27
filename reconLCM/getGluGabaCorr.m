
function GabaCorr = getGluGabaCorr(infile)

%%
%%
fileid = fopen(infile);
s = [];
while strcmp(s,'Correlation') == 0
    s = fscanf(fileid,'%s',1);
end
%%
s = fscanf(fileid,'%s',1);
%%
endTable = 0;
for index = 1:2
    temp = fscanf(fileid,'%s',1);
    disp(temp)
    corrMatCell{index} = temp;
end
%%
index = 3;
while endTable == 0
    temp = fscanf(fileid,'%s',1);
    if strcmp(temp,corrMatCell{1}) ~= 1
        ttemp = str2double(temp);
        if isnan(ttemp)
            corrMatCell{index} = temp;
        else
            corrMatCell{index} = ttemp;
        end
    else
        endTable = 1;
    end
    index = index + 1;
end
%%
corrMatCell = corrMatCell';
%%
GabaCorr = corrMatCell{93};
fclose(fileid);
