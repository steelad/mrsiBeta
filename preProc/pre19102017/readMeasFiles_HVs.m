%% Set up variables and read in files
function readMeasFiles_HVs(subjs,PATH)

for SUBJ = 1:length(subjs)
    cd(PATH)
    if isempty(subjs) == 1
        SNUM = 'subj'
    else
        SNUM = subjs{SUBJ};
        cd(SNUM);
    end
    %SNUM = ['HV_' num2str(subjs(SUBJ))]
    %FDIR = '/vols/Scratch/asteel/uzayData/';
    %cd([FDIR SNUM])
    FILES = dir('meas*slsr*');
    if size(FILES,1) == 0
        FILES = dir('meas*slaser*');
    end
    for i = 1:length(FILES)
	disp(['Detected: ' FILES(i).name])
    end

    %%
    outNameSTEM = [SNUM '_metab_water_dw_limitedvoi_']; %
    %addpath('/vols/Scratch/asteel/uzayScripts/');
    addpath('/vols/Scratch/asteel/uzayScripts/irt');
    %addpath('/vols/Scratch/asteel/uzayScripts/subject_1')
    setup;

    %%
    for FILE = 1:length(FILES);

    %%
    outName = [outNameSTEM num2str(FILE) '.mat'];

    %%
    curdir = pwd;
    FID = FILES(FILE).name;
    disp(['++ ' FID])
    FULLNAME = strcat(curdir,'/',FID);
    %%
    [data_water,refmrprot_water]=conceptd_data_read_kspace_24(FULLNAME);
    %%
    for ll=1:24
    for i=1:4
    temp_xx1_new (:,:,(ll-1)*4+i,:,1)=data_water(:,:,ll,:,i);
    end
    end
    %%
    for ll=1:24
    for i=5:8
        temp_xx1_new (:,:,(ll-1)*4+(i-4),:,2)=data_water(:,:,ll,:,i);
    end
    end
    %%

    nrings=24;
    FOV=240;
    ntheta  =   64;
    min_rad=0.25;
    ppr=64;
    k2  =  gen_k2_uzay(nrings,4,ntheta,FOV)*(FOV)/(nrings*2)*2*pi;
    k22=complex(k2(:,1),k2(:,2))';
    kk2_new = reshape(k22,[ppr,96]);
    %%
    Kx=real(kk2_new);
    Ky=imag(kk2_new);
    Kxx=repmat(Kx,[512,1]);
    Kyy=repmat(Ky,[512,1]);
    kx=pi/max(max(abs(Kxx)))*Kxx;
    ky=pi/max(max(abs(Kyy)))*Kyy;
    %%
    kf=linspace(-pi,pi,32768);  %% should be -pi to pi, just to normalize
    kf=kf';
    kkf=zeros(32768,96);
    for i=1:96
    kkf(:,i)=kf(:);
    end
    %%
    kf=kkf;  
    kx=[kx];
    ky=[ky];
    kf=[kf];
    kxx=kx(:);
    kyy=ky(:);
    kff=kf(:);
    %%
    st=nufft_init([kyy,kxx,kff], [24*2,24*2,512], [5,5,5], [4*24,4*24,512*2], [24,24,256],'minmax:kb');  % build the structure, remember to shift a little bit
    %%
    data=temp_xx1_new;
    for jjj=1:2
        disp(jjj)
    for j=1:32

    abc2=reshape(data(:,:,:,j,jjj),[32768,96]);

    data_a=[abc2];
    k_data=data_a; % apply density compensation
    Spect_data_2(:,:,:,j,jjj) = nufft_adj(k_data(:), st);


        end
    end
    %%

    for abc=1:2
        disp(abc)
     for j=1:32
        ll2=reshape(squeeze(Spect_data_2(:,:,:,j,abc)),[48*48,512]);
    ll3=fftshift(fft(ll2,[],2),2);
    ll4=reshape(ll3,[48,48,512]);
    fid_data(:,:,:,j,abc)=ll4;
          end
    end
    %load('uzay_water.mat','-mat');
    %%

    px = angle(squeeze((fid_data(:,:,1,:,2)+fid_data(:,:,1,:,1))/2));
    nws_water_nuf=squeeze((fid_data(:,:,1,:,2)+fid_data(:,:,1,:,1))/2).* exp( -1i * px );
    ll=estimate_csm_walsh(nws_water_nuf);
    ll_sq = sum(ll .* conj(ll),3); ll(ll < eps) = 1;
    %%
    for jj=1:512
    nws_water_nuf=squeeze((fid_data(:,:,jj,:,2)+fid_data(:,:,jj,:,1))/2).* exp( -1i * px );
    img_water(:,:,jj)= sum(conj(ll) .*nws_water_nuf ,3) ./ ll_sq;
    end
    %%


    for kk=1:2
        disp(kk)
    for jj=1:512
          nws_metab_nuf=squeeze(fid_data(:,:,jj,:,kk)).* exp( -1i * px );
      img_metab(:,:,jj,kk) = sum(conj(ll) .* nws_metab_nuf,3) ./ ll_sq;
     end
    end
    %%

    save(outName,'img_metab','refmrprot_water','img_water','FID','-mat');
    end
end
