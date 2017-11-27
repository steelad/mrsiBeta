function segT1Func(subj,PATH,filename)
if  ~isempty(PATH) && ~strcmp(PATH,'.')
    cd([PATH,subj,'/T1'])
end
setenv('f',filename)
!echo $f
tic
!echo Stripping $f
!3dSkullStrip -input $f -prefix T1_ss.nii
!echo Segmenting $f
!fast -p T1_ss.nii
toc
