function visualize_LCModel(subj,filename,PATH);
%%

%%
%subj = 'HV_105';
%%';
if isempty(filename) ~=0
    load(['../' subj,'_finalWksp512_noAvg.mat'],'dw_mean_hlsvd','meta_mask')
else
    load(filename,'dw_mean_hlsvd','meta_mask')%,'postLipid_img'%% FIX THIS WHEN COMPLETE 18012018
end
%%
if isempty(PATH)==0
    disp(PATH)
    cd(PATH)
end
%% get file size params
inform = size(dw_mean_hlsvd);

if size(inform,2) < 4
    warning('Only single timepoint detected.')
    inform(4) = 1;
end
% try
%     isempty(inform(4))
% catch
%     inform(4) = 1;
% end
%% determine subject name length( relative to delimiter)
temp1 = length(strsplit(subj,'_'));
tempf = dir('*_lcm');
tempf = strsplit(tempf(1).name,'_');
%if inform(4) == 1
%%    stem = strjoin(tempf(temp1+3:end),'_');
%else
    stem = strjoin(tempf(temp1+4:end),'_');
%end
%%
stem
% %%
 tempMask = meta_mask;
% tempMask(29,32) = 0;
%%
tt = [];
for i = size(tempMask,2):-1:1;
    temp = tempMask(:,i);
    tt = find(temp == 1);
    if isempty(tt) == 0
        break
    end
 end
% %% adjust mask as needed
% tempMask(min(tt),:) = zeros(1,size(tempMask,2)); 
% tempMask(max(tt),:) = zeros(1,size(tempMask,2)); 
% tempMask(:,[i-1 1]) = zeros(size(tempMask,1),2);

%%
CC= 1;
 for i=1:inform(1)
     for k=1:inform(2)
         if meta_mask(i,k) ~=0
            for MM = 1:inform(4)
                 
                B=[ subj '_' int2str(i) '_' int2str(k) '_' num2str(MM) '__' stem '/lcm_1.COORD'];
                if ~exist(B)
                    B = [ subj '_' int2str(i) '_' int2str(k) '_' stem '/lcm_1.COORD'];
			%pwd
			%B
                    assert(exist(B,'file')== 2,'++ File not present. Something is wrong')
                end
               
            try
             [lcmodelresults,checkFATAL]=readcoord(B);

                 for ii=1:length(lcmodelresults.metabconc)
                     if  checkFATAL ~= 1
                         if CC == 1
                             
                            metab_map=zeros(48,48,length(lcmodelresults.metabconc)+2,inform(4));
                            metab_map2=zeros(48,48,length(lcmodelresults.metabconc)+2,inform(4));
                            metab_map1=zeros(48,48,length(lcmodelresults.metabconc),inform(4));
                            gluGabaCorrMap=zeros(48,48,inform(4));
                            lcm_output_struct = cell(48,48,inform(4));
                            CC = 0;
                            
                         end
                     lcm_output_struct{i,k,MM} = lcmodelresults;    
                     metab_map(i,k,ii,MM)=lcmodelresults.metabconc(ii).absconc;
                     metab_map2(i,k,ii,MM)=lcmodelresults.metabconc(ii).relconc;
                     metab_map1(i,k,ii,MM)=lcmodelresults.metabconc(ii).SD;
                     else
                          metab_map(i,k,:,MM)=0;
                     end
                 end
                if checkFATAL ~= 1
                    C=[ subj '_' int2str(i) '_' int2str(k) '_' num2str(MM) '__' stem '/lcm_1.PRINT'];
                if ~exist(C)
                    C = [ subj '_' int2str(i) '_' int2str(k) '_' stem '/lcm_1.PRINT'];
                    assert(exist(B,'file')== 2,'++ File not present. Something is wrong')
                end
                    corrGabaGlu = getGluGabaCorr(C);
                    gluGabaCorrMap(i,k,MM) = corrGabaGlu;

                    metab_map(i,k,length(lcmodelresults.metabconc)+1,MM)=lcmodelresults.SN;
                    metab_map(i,k,length(lcmodelresults.metabconc)+2,MM)=lcmodelresults.linewidth;
                    metab_map2(i,k,length(lcmodelresults.metabconc)+1,MM)=lcmodelresults.SN;
                    metab_map2(i,k,length(lcmodelresults.metabconc)+2,MM)=lcmodelresults.linewidth;
                end


             catch
                 disp(['++ ' B ' Fail.'])
             end
             end

         end
     end
end

 %%
 
%%
%%
% figure;
% imagesc((real(gluGabaCorrMap(:,:,1))))
% axis square
% axis off
% 
% %%
% figure;
% imagesc((real((metab_map2(:,:,8,1)))))
% axis square
% axis off
% %%
% figure;
% imagesc(fliplr(real(imrotate(metab_map2(:,:,5),116))))
% axis square
% axis off
%%

save([subj '_metab_maps.mat'],'metab_map','metab_map1','metab_map2','gluGabaCorrMap','lcm_output_struct')
%%
