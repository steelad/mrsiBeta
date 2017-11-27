function bc = issameclass(classtype)

%It gives a matrix bc of size length(classtype) x length(classtype)
%bc(i,j)=1 if classtype of voxel i = classtype of voxel j

lc = length(classtype);
bc=zeros(lc,lc);

% checking classe types
if isa(classtype,'double') % the classes are represented by numbers
    tc = 1;
elseif iscell(classtype) % the classes are represented by strings
    tc = 0;
else
    errordlg('The variables has to be a vector of integers (format double) or a cell of strings');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tc==1 % the classes are represented by numbers
for i=1:lc
    for j=1:lc
        if classtype(i)==classtype(j)
            bc(i,j)=1;
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tc==0 % the classes are represented by strings
for i=1:lc
    for j=1:lc
        if strcmp(classtype(i),classtype(j))
            bc(i,j)=1;
        end
    end
end
end