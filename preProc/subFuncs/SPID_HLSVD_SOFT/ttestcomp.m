function ttestcomp(classresfilename) 
% t-test for each metabolite of each pair of classes
v = load(classresfilename);
[classopt,procres] = checkclassvar(a);

classtype = classopt.classtype;
scores = procres.scores;

% list of class labels 
clabels = classlabels(classopt);
lcl = length(clabels);
for i=1:lcl
    class(i).matrix = [];
end

%creating class matrices
for i = 1:size(scores,1)
    if iscell(clabels)
        ind = find(strcmp(clabels,classtype(i)));
    else
        ind = find(clabels==classtype(i));       
    end
    class(ind).matrix = [class(ind);scores(i,:)]; 
end

comb = factorial(lcl)/(factorial(2)*factorial(lcl-2));
c = 1;
pc = [];
for i=1:comb
    j=i;
    while j<=comb
        pc = [pc;i j];
        j=j+1;
    end
end
    
alpha = 0.05;
% testing all the null hypothesis 
for i = 1:comb % i possible pair of classes 
    for j=1:size(scores,2)   %j features (ex: metabolites)
        h(j,i) = ttest2(class(pc(i,1)).scores(:,j),class(pc(i,2)).scores(:,j),alpha); 
    end
end

%displaying the results 
disp('T-test comparisons: H0: mean(metabolite_i of class j1)==mean(metabolite_i of class j2)')
for i = 1:comb
    disp(['classes ',clabels(pc(i,1)), 'and ', clabels(pc(i,2))])
    disp(str2num(mean(class(pc(i,1)))'), '+- ',str2num(std(class(pc(i,1)))'), '    ', ...
        str2num(mean(class(pc(i,2)))'), '+- ',str2num(std(class(pc(i,2)))'), '     ', h(:,i));
end
