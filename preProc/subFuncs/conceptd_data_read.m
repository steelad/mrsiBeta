
function [data,refmrprot] = concept_data_read(filename)
[refmrprot, refmdh, inref] =rdMeas_dene(filename);
ppr = 64; % pts per ring (oversampled)
%data = reshape(data,ppr,t2pts,phases,channels);
inref_n=inref(9:end,:,:);
t2pts=(size(inref_n,1))/ppr;
min_rad=0.25;
gd1 = 1.5; % def is 1.5
gd2 = 1; % def is 1

rep=refmrprot.sWiPMemBlock.alFree (2);

for i=1:t2pts
    temp(:,i,:,:)=inref_n((i-1)*ppr+1:i*ppr,:,:);
   
      temp2(:,i,:,:)=circshift(inref_n((i-1)*ppr+1:i*ppr,:,:),8,1);
   
    %data(:,i,:,:) = inref_n((i-1)*ppr+1:i*ppr,:,:);
end




for ll=1:16
for jj=1:20
for i=1:32
    temp_xx1(:,:,jj,i,ll)=[temp(:,:,i+(jj-1)*64+1280*(ll-1)) circshift(circshift(temp(:,:,i+32+(jj-1)*64+1280*(ll-1)),8,1),0,1)];
end
   
end
end
% 
% for ll=2:2:10
% 
% 
%     temp_xx1(:,:,:,:,2)=circshift(temp_xx1(:,:,:,:,ll),32,1);;
% 
%  
% end

 data=temp_xx1;
 
    
    return;