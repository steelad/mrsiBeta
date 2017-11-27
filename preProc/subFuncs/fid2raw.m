
function complexfid=fid2raw(filename,vivofid,index)





fid(:,1)=real(vivofid);
fid(:,2)=imag(vivofid);
directory=pwd;

newdir=[filename '_lcm']
mkdir(newdir);

fileid=fopen([newdir '/lcm_' num2str(index) '.RAW'],'w');
            fprintf(fileid,[' $NMID ID=''' 'lcm_' num2str(index) ''', FMTDAT=''(8E13.5)''\n']);
            fprintf(fileid,[' TRAMP=1., VOLUME=1. $END\n']);
            fprintf(fileid,'%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n',[fid(:,1)'; fid(:,2)']); 
            fclose(fileid);

            
end