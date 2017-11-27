%PC plotting

cdir = getappdata(0,'currentdir');
[fn p]=uigetfile({'*.mat','MAT-files (*.mat)';},'Load file from');
if isequal(fn,0) | isequal(p,0)
    return
end
%cd(pwd);
prompt={'How many components do you want to plot (2 or 3)'};
name='Input for PC plot';
numlines=1;
defaultanswer={'2'};
answer=inputdlg(prompt,name,numlines,defaultanswer);
if isempty(answer)
    disp('...cancelled')
    return
end
npcs = [str2num(char(answer(1)))];

a = load(fullfile(p,fn));
[classopt,procres] = checkclassvar(a);

[k,comb] = binaryencoding(classopt.classtype);
clabels = classlabels(classopt);

for i =1:comb%comb
        ind = find(k(:,i) ~= 0);
        if length(size(procres.scores))==3
            data = procres.scores(:,:,i);
        else
            data = procres.scores;
        end
        ind1 = find(k(:,i) == -1);
        ind2 = find(k(:,i) == 1);
        [pc, zscores, pcvars] = princomp(data(ind,:),'econ');
        figure(100+i)
        if npcs==3
        plot3(zscores(ind1,1),zscores(ind1,2),zscores(ind1,3),'b*')
        hold on 
        plot3(zscores(ind2,1),zscores(ind2,2),zscores(ind2,3),'r+')
        else
        plot(zscores(ind1,1),zscores(ind1,2),'b*')
        hold on 
        plot(zscores(ind2,1),zscores(ind2,2),'r+')
        end
        grid on
end
