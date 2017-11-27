function [signal, ndp, step, nfids, volume]= read_philips_file(file2read)

fin	= fopen([(file2read) '.spar'],'rt');
tmsrf	= 0;
while 1
    [dummy,count] = fscanf(fin,'%s',1);
    if count < 1
        break;
    end
    if strcmp(dummy,'nucleus') == 1
        dummy = fscanf(fin,'%s',1);
        dummy = fscanf(fin,'%s',1);
        if strcmp(dummy,'"1H"') | strcmp(dummy,'"1h"') | strcmp(dummy,'1H')
            nuc = 'h';
        elseif strcmp(dummy,'"31P"') | strcmp(dummy,'"31p"') | strcmp(dummy,'31P')
            nuc = 'p';
        elseif strcmp(dummy,'"13C"') | strcmp(dummy,'"13c"') | strcmp(dummy,'13C')
            nuc = 'c';
        elseif strcmp(dummy,'"19F"') | strcmp(dummy,'"19f"') | strcmp(dummy,'19F')
            nuc = 'f';
        else
            nuc = dummy(2:length(dummy)-1);
        end
    elseif strcmp(dummy,'synthesizer_frequency') == 1
        dummy    = fscanf(fin,'%s',1);
        sfrq	= fscanf(fin,'%f',1);
    elseif strcmp(dummy,'offset_frequency') == 1
        dummy    = fscanf(fin,'%s',1);
        woff	= fscanf(fin,'%f',1);
    elseif strcmp(dummy,'sample_frequency') == 1
        dummy    = fscanf(fin,'%s',1);
        sw	= fscanf(fin,'%f',1);
    elseif strcmp(dummy,'samples') == 1
        dummy    = fscanf(fin,'%s',1);
        ndp	= fscanf(fin,'%d',1);
    elseif strcmp(dummy,'averages') == 1
        dummy    = fscanf(fin,'%s',1);
        scans	= fscanf(fin,'%d',1);
    elseif strcmp(dummy,'ap_size') == 1
        dummy = fscanf(fin, '%s',1);
        ap_size= fscanf(fin, '%d',1);
    elseif strcmp(dummy,'lr_size') == 1
        dummy = fscanf(fin, '%s',1);
        lr_size= fscanf(fin, '%d',1);
    elseif strcmp(dummy,'cc_size') == 1
        dummy = fscanf(fin, '%s',1);
        cc_size= fscanf(fin, '%d',1);
    elseif strcmp(dummy,'nr_time_series_spectra') == 1
        dummy    = fscanf(fin,'%s',1);
        nfids	= fscanf(fin,'%d',1);
        if nfids < 1; nfids	= 1; end
        tmsrf	= 1;
    end
    if tmsrf == 0
        if strcmp(dummy,'spec_num_row') == 1
            dummy    = fscanf(fin,'%s',1);
            nfids	= fscanf(fin,'%d',1);
        end
    end
end
fclose(fin);


%Rene zegt: bij philips begin keihard op 0 zetten!

begin=0;


%gedeeld door 1000 omdat in de .spar file de samplefreq gegeven is in Hz, terwijl
%MRUI en zo met kHz werkt.
step	= 1000.0 / sw;

pmsdat=[(file2read) '.sdat'];
pmspar=[(file2read) '.spar'];

bpar=[1000 1 0];


pmsfid       = zeros(nfids,ndp);
fin = fopen(pmsdat,'r','d');
for scnt = 1:nfids
    [sig1,s1cnt]       = fread(fin,[2 ndp],'float32');
    pmsfid(scnt,1:ndp) = bpar(2) * (sig1(1,1:ndp) + i *sig1(2,1:ndp));
end
fclose(fin);
ftit         = 'Philips data file, converted format (*.sdat)';

signal=pmsfid;

volume=ap_size*lr_size*cc_size/1000;  %volume in ml!!!!!
