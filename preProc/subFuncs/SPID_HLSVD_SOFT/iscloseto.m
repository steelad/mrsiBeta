function bs = iscloseto(position,opt)

% this function generates a matrix bs of size LxL where L is
% size(position,1), i.e. the number of signals or voxels
% the elements of the matrix are boolean, s.t. bs(i,j)=1 if voxel i is
% close to voxel j and bs(i,j)=0 otherwise
% opt=0 => only the adjacent voxels (4 voxels:top, bottom, left and right) are considered as close
% opt=0 => the adjacent voxels + voxels on the diagonal are considered as
% close (=8 voxels)
sp = size(position,1);
mr = max(position(:,1));
mc = max(position(:,2));
bs = zeros(sp,sp);
for i=1:sp
    al = [position(i,1)-1,position(i,2);position(i,1),position(i,2)+1;position(i,1)+1,position(i,2);position(i,1),position(i,2)-1;...
        position(i,1)-1,position(i,2)-1;position(i,1)-1,position(i,2)+1;position(i,1)+1,position(i,2)+1;position(i,1)+1,position(i,2)-1];
    inddel = [];
    if position(i,1)==1
        inddel = [inddel 1 5 6];
    end
    if position(i,1)==mr
        inddel = [inddel 3 7 8];
    end
    if position(i,2)==1
        inddel = [inddel 4 5 8];
    end
    if position(i,2)==mc
        inddel = [inddel 2 6 7];
    end
    al(inddel,:) = zeros(length(inddel),2);
    ind=[];
    if opt
        for j = 1:8
            ind = [ind (al(j,2)-1)*mr+al(j,1)];
        end
        ind(ind==0)=[];
    else
        al = al(1:4,:);
        for j = 1:4
            ind = [ind (al(j,2)-1)*mr+al(j,1)];
        end
        ind(ind<=0)=[];
    end
    bs(i,ind)=ones(1,length(ind));
end
