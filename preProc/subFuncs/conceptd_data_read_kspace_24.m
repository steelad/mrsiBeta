
function [data,refmrprot] = concept_data_read_kspace(filename)
[refmrprot, refmdh, inref] =rdMeas_dene(filename);
ppr = 64; % pts per ring (oversampled)
inref_n=inref(9:end,:,:);
t2pts = (length(inref_n))/ppr;
% min_rad=0.25;
% gd1 =5; % def is 1.5
% gd2 = 2; % def is 1

%rep=refmrprot.sWiPMemBlock.alFree (2);
for i=1:t2pts
    temp(:,i,:,:)=inref_n((i-1)*ppr+1:i*ppr,1:2:end,:);
      temp2(:,i,:,:)=circshift(inref_n((i-1)*ppr+1:i*ppr,2:2:end,:),8,1);
   
    %data(:,i,:,:) = inref_n((i-1)*ppr+1:i*ppr,:,:);
end

for i=1:8
    temp_new(:,:,:,:,i)=temp(:,:,(i-1)*24+1:i*24,:);
     temp_new2(:,:,:,:,i)=temp2(:,:,(i-1)*24+1:i*24,:);
end
 data=[temp_new temp_new2];
 
return;