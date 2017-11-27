function complexfid=runlcmodel_ch_7T(filename,TEB)


load([filename '.spa'],'-mat')
%dic_name=[dirs(3).name(abc_dcm(end)+1:end)]
nbpoints=length(data.metab); % Varian np divided by 2
ntmetab = data.ntmetab;
%load('ch.mat')
sw =  data.params(2);
sw_hz=sw;
sfrq = data.params(1);
TE=14;
 T2metab = 100;
 T2water = 50;  %=(70*1.15)
tissuewatercontent = 1;  % tissue water content may vary depending on the brain region. (eg 0.71 for pure white matter, 0.81 for pure grey matter)
CSFwatercontent = 1;      % assuming CSF is pure water
attmet = 1.5;
	% number of averages in water ref file
concpurewater = 55555;  % concentration of pure water in mM
%data.ntws=1;
atth2o = 1;
%if (waterTEs.CSFfrac<0)
    CSFfrac=0;
%else
 %   CSFfrac=waterTEs.CSFfrac/100;
%end

            %factor of 1.15 from Michaeli MRM 2003 estimated at 3T
%atth2o = data.ntws;%comment
wconc = concpurewater*(CSFfrac*CSFwatercontent + (1-CSFfrac)*tissuewatercontent);
    wconc = wconc./(1-CSFfrac);

       fid2raw(filename,data.metab,1);
       fid2h2o(filename,data.ws,1);
       
       
       
       
    

    tes_basis(1).te = '/home/fs0/asteel/scratch/uzayScripts/scripts_templates/scripts/mrsiBeta/misc/3T_slaser_32vespa_1250';

%7T_real_newmac_2013_oxford_essen_28ut.BASIS
%tes_basis(1).te='/home/orochi3-raid7/uzay-data/mc7t/Magnetom7T-SVS_SLASER30T7+2HG+MM_CMRR'
%/home/rosalind-raid1/uzay/Metab_Files/7T/2012_11_15_essen_22ut
%/home/rosalind-raid1/uzay/Metab_Files/7T/2012_11_15_essen_19ut
%/home/rosalind-raid1/uzay/Metab_Files/7T/2012_11_16_essen_25ut/7T_sead_mac_2012_nih_essen_25ut.BASIS
%/home/rosalind-raid1/uzay/Metab_Files/7T/2012_12_20_mc7t_cmrr_28ut/7T_sead_mac_2012_cmrr_essen_28ut.BASIS
%/home/rosalind-raid1/uzay/Metab_Files/7T/2013_01_09_philips_298/7T_sead_mac_2013_philips_essen_28ut

directory=pwd;

newdir=[filename '_lcm']
mkdir(directory,newdir);
lcmodeloutputfilename=['lcm_' num2str(1)];
%fileid=fopen([lcmodeloutputfilename '.CONTROL'],'w');

fileid=fopen([newdir '/lcm_' num2str(1) '.CONTROL'],'w');
fprintf(fileid,' $LCMODL\n');
fprintf(fileid,[' TITLE=''' lcmodeloutputfilename '''\n']);
fprintf(fileid,[' OWNER=''FMRIB Centre, University of Oxford''\n']);
fprintf(fileid,' KEY=217638264\n');  % Key for Bernoulli
fprintf(fileid,[' PGNORM=''US''\n']);
fprintf(fileid,[' FILPS=''' lcmodeloutputfilename '.PS''\n']);
fprintf(fileid,[' FILCOO=''' lcmodeloutputfilename '.COORD''\n']);


fprintf(fileid,' NCOMBI=5\n');
fprintf(fileid,[' CHCOMB(1)=''PCho+GPC''\n']);
fprintf(fileid,[' CHCOMB(2)=''Cr+PCr''\n']);
fprintf(fileid,[' CHCOMB(3)=''NAA+NAAG''\n']);
fprintf(fileid,[' CHCOMB(4)=''Glu+Gln''\n']);
fprintf(fileid,[' CHCOMB(5)=''Glc+Tau''\n']);
fprintf(fileid,[' NAMREL=''Cr+PCr''\n']);
fprintf(fileid,' CONREL=8.00\n');
fprintf(fileid,' LPRINT=6\n');
fprintf(fileid,[' FILPRI=''' lcmodeloutputfilename '.PRINT''\n']);
fprintf(fileid,' LCOORD=9\n');
fprintf(fileid,' NOMIT=17\n'); %for more CHOMIT change NOMIT from 1 to 7
fprintf(fileid,[' CHOMIT(1)=''Ace''\n']);
fprintf(fileid,[' CHOMIT(2)=''ATP''\n']);
fprintf(fileid,[' CHOMIT(3)=''Gly''\n']);
fprintf(fileid,[' CHOMIT(4)=''Thr''\n']);
fprintf(fileid,[' CHOMIT(5)=''Cho''\n']);
fprintf(fileid,[' CHOMIT(6)=''Gua''\n']);
fprintf(fileid,[' CHOMIT(7)=''Ser''\n']);
fprintf(fileid,[' CHOMIT(8)=''-CrCH2''\n']);
fprintf(fileid,[' CHOMIT(9)=''Lip13a''\n']);
fprintf(fileid,[' CHOMIT(10)=''Lip13b''\n']);
fprintf(fileid,[' CHOMIT(11)=''Lip13c''\n']);
fprintf(fileid,[' CHOMIT(12)=''Lip13d''\n']);
fprintf(fileid,[' CHOMIT(13)=''Lip09''\n']);
%fprintf(fileid,[' CHOMIT(14)=''Lip20''\n']);
fprintf(fileid,[' CHOMIT(15)=''Mac''\n']);
fprintf(fileid,[' CHOMIT(16)=''bHB''\n']);
fprintf(fileid,[' CHOMIT(17)=''HG''\n']);
fprintf(fileid,' NUSE1=5\n');
fprintf(fileid,[' CHUSE1(1)=''Cr''\n']);
fprintf(fileid,[' CHUSE1(2)=''PCr''\n']);
fprintf(fileid,[' CHUSE1(3)=''NAA''\n']);
fprintf(fileid,[' CHUSE1(4)=''Glu''\n']);
fprintf(fileid,[' CHUSE1(5)=''Ins''\n']);
fprintf(fileid,' NCALIB=0\n');
fprintf(fileid,' WSPPM=3.0241\n');
fprintf(fileid,' WSMET=''Cr''\n');



 
 fprintf(fileid,' NSIMUL=14\n');
fprintf(fileid,' NRATIO=0\n');
 fprintf(fileid,' NEACH=99\n');
 

      




fprintf(fileid,[' DEGZER=0.0\n']);
fprintf(fileid,' SDDEGZ=5.0\n');
fprintf(fileid,[' DEGPPM=0.\n']);
fprintf(fileid,[' SDDEGP=5.\n']);
fprintf(fileid,' SHIFMN=-0.2,-0.1\n');
fprintf(fileid,' SHIFMX=0.3,0.3\n');
fprintf(fileid,' FWHMBA=0.005\n');
fprintf(fileid,' RFWHM=2.50\n');
fprintf(fileid,' DKNTMN=0.25\n');
fprintf(fileid,' PPMST=4.2\n');
fprintf(fileid,' PPMEND=0.5\n');

%%%
fprintf(fileid,' DOWS=T\n');
fprintf(fileid,' DOECC=F\n');
fprintf(fileid,[' ATTH2O=' num2str(atth2o*1) '\n']);

fprintf(fileid,[' ATTMET=' num2str(attmet) '\n']);

fprintf(fileid,[' WCONC=' num2str(wconc) '\n']);
%%%%
%fprintf(fileid,' NKEEP=0\n');
%fprintf(fileid,' CHKEEP=0\n');

%sfrq=123.1788
%%%%%%
fprintf(fileid,[' FILRAW=''' lcmodeloutputfilename '.RAW''\n']);
fprintf(fileid,[' FILH2O=''' lcmodeloutputfilename '.H2O''\n']);
%fprintf(fileid,[' FILBAS=''4T_SIM_dinesh10ms.BASIS''\n']);
fprintf(fileid,[' FILBAS=''' tes_basis(0+1).te '.BASIS''\n']);
fprintf(fileid,[' DELTAT=' num2str(1/sw_hz) '\n']);
fprintf(fileid,[' NUNFIL=' num2str(nbpoints) '\n']);%%nbpoints should be defined
fprintf(fileid,[' HZPPPM=' num2str(sfrq) '\n']);
%fprintf(fileid,' VITRO=T\n');
fprintf(fileid,' PPMSHF=0.0\n');

fprintf(fileid,' $END\n');
fclose(fileid);
folders=[newdir]

%unix(['ssh lcmodel " cd '  folders ';/home/fs0/uemir/.lcmodel/bin/lcmodel < ' lcmodeloutputfilename '.CONTROL"'],'-echo')

cd(folders)
unix(['/home/fs0/uemir/.lcmodel/bin/lcmodel < ' lcmodeloutputfilename '.CONTROL'])
cd(directory)
