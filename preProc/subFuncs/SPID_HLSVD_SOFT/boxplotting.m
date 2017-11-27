%boxplotting

cdir = getappdata(0,'currentdir');
[fn p]=uigetfile({'*.mat','MAT-files (*.mat)';},'Load file from');
if isequal(fn,0) | isequal(p,0)
    return
end
%cd(pwd);

a = load(fullfile(p,fn));
[classopt,procres] = checkclassvar(a);

[k,comb] = binaryencoding(classopt.classtype);
clabels = classlabels(classopt);
c1 = 1;
c2 = 2;
[sp1,sp2] = size(procres.scores);
v = version;
if str2num(v(1:3))<7.1
    errordlg('Matlab version 7.1 and more recent is needed to use this tool')
    return
end

for i =1:comb%comb
    ind = find(k(:,i) ~= 0);
    if length([sp1,sp2])==3
        data = procres.scores(:,:,i);
    else
        data = procres.scores;
    end
    ind1 = find(k(:,i) == -1);
    ind2 = find(k(:,i) == 1);
    if isfield(procres,'misc')
        if isfield(procres.misc,'compNames')
            cl = procres.misc.compNames;
        else
            for j = 1:sp2
                cl(j) = {num2str(j)};
            end
        end
    else
        for j = 1:sp2
            cl(j) = {num2str(j)};
        end
    end
    figure(100+i)
    boxplot(data(ind2,:),'position',1.2:1:sp2-1+1.2,'widths',0.1,'colors','g');
    hold on;
    boxplot(data(ind1,:),'position',0.8:1:sp2-1+0.8,'widths',0.1,'labels',cl);
    if isa(clabels(c1),'cell')
        title(['Boxplot class ', char(clabels(c1)), ' vs ', char(clabels(c2))])
    else
        title(['Boxplot class ', num2str(clabels(c1)), '(blue) vs ', num2str(clabels(c2)), ' (green)'])
    end
    %axis(0:sp2+1)
    if c2 == length(clabels)
        c1 = c1+1;
        c2 = c1+1;
    else
        c2 = c2+1;
    end
end