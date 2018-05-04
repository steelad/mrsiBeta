% read .COORD file and plot data
% PGH Oct 2001
%
% verify version of LCModel used  - DKD April 2012
%
% INPUT PARAMETER : filename
% OUTPUT PARAMETER : lcmodelresults structure
%       spectrumppm             : ppm values
%       spectrumdata            : raw data before fit
%       spectrumfit             : fit
%       spectrumbasl            : baseline
%       metabconc               : structure containing output information on fitted concentrations
%               name            : name of metabolite
%               relconc         : relative concentration
%               absconc         : absolute concentration
%               SD              : standard deviation (Cramer-Rao Bound)
%       linewidth               : estimated linewidth
%       SN                      : estimated signal-to-noise ratio

function [lcmodelresults,checkFATAL]=readcoord(filename)

% define output structure for results

lcmodelresults=struct('spectrumppm',0,'spectrumdata',0,'spectrumfit',0,'spectrumbasl',0,'metabconc',0,'linewidth',0,'SN',0);
lcmodelresults.metabconc=struct('name',0,'relconc',0,'absconc',0,'SD',0);

% open .COORD file

fprintf('**** READING LCMODEL RESULTS IN .COORD FILE *****\n\n');
fprintf(['Opening file ' filename '\n']);
fileid=fopen(filename);
%return if .COORD does not exist
if (fileid == -1)
    fprintf(['Cannot find ' filename '\n EXITING\n']);
    return
end
checkFATAL = 0;
s = [];
while strcmp(s,'input') == 0
    s = fscanf(fileid,'%s',1);
    if strcmp(s,'FATAL') == 1
        checkFATAL = 1;
    end
end
fclose(fileid);
if checkFATAL == 1
    disp(['++ Fatal error detected.'])
    return
end
fileid=fopen(filename);
% discard text until "Version" is found
s=[];
while strcmp(s,'(Version')==0
    s=fscanf(fileid,'%s',1);
end
s=fscanf(fileid,'%s',1);
versionID = str2double(s(1:3));

% discard text until beginning of concentrations table (preceded by word 'Metab.')
s=[];
while strcmp(s,'Metabolite')==0
    s=fscanf(fileid,'%s',1);
end

% read concentration values
index=1;
endtable=0;
while (endtable==0)
    lcmodelresults.metabconc(index).absconc=fscanf(fileid,'%f',1);
    lcmodelresults.metabconc(index).SD=fscanf(fileid,'%f',1);
    temp=fscanf(fileid,'%s',1);                  % read and discard '%' character
    if strcmp(temp,'lines'), endtable=1; end   % if word 'lines' found then concentration table has been completely read
    lcmodelresults.metabconc(index).relconc=fscanf(fileid,'%f',1);
    lcmodelresults.metabconc(index).name=fscanf(fileid,'%s',1);
    index=index+1;
end

lcmodelresults.metabconc=lcmodelresults.metabconc(1:length(lcmodelresults.metabconc)-1); % discard last line of table
fprintf([num2str(length(lcmodelresults.metabconc)) ' metabolite concentrations values have been read\n'])

% discard text until linewidth (preceded by word 'FWHM')

s=[];
while strcmp(s,'FWHM')==0
    s=fscanf(fileid,'%s',1);
end

% read linewidth

fscanf(fileid,'%s',1); % discard '='
lcmodelresults.linewidth=fscanf(fileid,'%f',1);

% discard text until S/N (preceded by word 'S/N')

s=[];
while strcmp(s,'S/N')==0
    s=fscanf(fileid,'%s',1);
end

% read S/N

fscanf(fileid,'%s',1); % discard '='
lcmodelresults.SN=fscanf(fileid,'%f',1);

% discard text until number of data points (preceded by word 'extrema')

s=[];
while isempty(findstr(s,'extr'))
    s=fscanf(fileid,'%s',1);
end

% read number of points

nbpoints=fscanf(fileid,'%d',1);

% read and discard text 'points on ppm-axis = NY'

fscanf(fileid,'%s',5);

% read ppm values

lcmodelresults.spectrumppm=fscanf(fileid,'%f',nbpoints);
fprintf([num2str(length(lcmodelresults.spectrumppm)) ' ppm values have been read\n'])

% read and discard text 'NY phased data points follow'

fscanf(fileid,'%s',5);

% read data values

lcmodelresults.spectrumdata=fscanf(fileid,'%f',nbpoints);
fprintf([num2str(length(lcmodelresults.spectrumdata)) ' data values have been read\n'])

% read and discard text 'NY points of the fit to the follow'
if (versionID>=6.2)
    fscanf(fileid,'%s',9);
else
    fscanf(fileid,'%s',8);
end

% read fit values

lcmodelresults.spectrumfit=fscanf(fileid,'%f',nbpoints);
fprintf([num2str(length(lcmodelresults.spectrumfit)) ' fit values have been read\n'])

% read and discard text 'NY background values follow'

fscanf(fileid,'%s',4);

% read baseline values

lcmodelresults.spectrumbasl=fscanf(fileid,'%f',nbpoints);
fprintf([num2str(length(lcmodelresults.spectrumbasl)) ' baseline values have been read\n\n'])

% close .COORD filec
% this is super hacky.
lcmodelresults.metFit = struct;
c = 0;
%for i = 1:length(lcmodelresults.metabconc)
cmplt = '';
metLen = length(lcmodelresults.spectrumbasl);
while metLen ==  length(lcmodelresults.spectrumbasl)
    c = c+1;
    met = fscanf(fileid,'%s',1);
    disp(met)
    fscanf(fileid,'%s',3); %discard unnecessary bits%
    
    tempMet = fscanf(fileid,'%f',nbpoints);
    metLen = length(tempMet);
    
    lcmodelresults.metFit(c).name = met;
    lcmodelresults.metFit(c).modFit = tempMet ;
    
    fprintf([num2str(length(lcmodelresults.metFit(c).modFit)) ' ' ...
       lcmodelresults.metFit(c).name ' values have been read\n\n'])
end
lcmodelresults.metFit = lcmodelresults.metFit(1:end-1) ;% clean up.
fclose(fileid);

