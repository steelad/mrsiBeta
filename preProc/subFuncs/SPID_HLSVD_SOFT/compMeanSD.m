%***************************************************************************
%                                compMeanSD.M
%***************************************************************************
% PURPOSE:  computing the mean and the standard deviation for a paricular
%           tumour
%***************************************************************************
% CALL:     compMeanSD(protocol,tumour,resultf)
%***************************************************************************
% INPUT:    tumour      =  name of the tumour
%           protocol    = number of the protocol 
%           resultsf    = result file
% OUTPUT:   mean        = mean of the set of data (vector, size = number of metabolites*4)
%                           4 = number of parameters per metabolites (amplitude,phase,damping,frequency)                    
%***************************************************************************
% COMMENTS: resultf must contains the parameters: ampl,damp,freq,phase,
%           compNames
%           where 
%           ampl,damp,freq,phase are matrix vectors (rows=signals, columns=metabolites)
%           compNames contains the names of the metabolites (column vector)
%***************************************************************************
function compMeanSD(protocol,tumour,resultf) 
global ampl damp freq phase
global compNames


%OS 
if isunix
    OS='/';
else
    OS='\';
end

lCN = length(compNames);

%LOADING THE RESULT FILES
eval(['load ',resultf])

% row vectors for means and standard deviations 
meanampl = mean(ampl);
meandamp = mean(damp);
meanfreq = mean(freq);
meanphas = mean(phase);

SDampl = std(ampl);
SDdamp = std(damp);
SDfreq = std(freq);
SDphas = std(phase);

% One row per metabolite (first 4 columns = means, last 4 c. = SD)
means = [meanampl' meanphas' meandamp' meanfreq'];
stds = [SDampl' SDphas' SDdamp' SDfreq'];
meanstd = [means stds];
meanstd = meanstd';

% file to save to
savetof = ['..',OS,protocol,OS,'tables',OS,tumour,'.txt'];
compNames = cellstr(compNames)

fid = fopen(savetof,'w');
fidcode='\n';
for i=1:lCN
fidcode = [fidcode,' %f \t\t'];
end
fidcode2 = ['\n %f \t\t %f \t\t %f \t\t %f \t\t %f \t\t %f \t\t %f \t\t %f'];
for i=1:lCN
fprintf(fid,'%s \t\t',char(compNames(i)));
end
fprintf(fid,fidcode2,meanstd);
fclose(fid);
