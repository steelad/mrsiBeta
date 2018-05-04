function reconOut = hlsvdProc_SVMRS(img1,timeStep,prefix)
%%

disp('++ HLSVD being applied to all images in set... This will take a while.')
%%
disp(size(img1))
numImgs = size(img1,2);
tps = size(img1,1);
tic
reconOut = zeros(size(numImgs));

for ii=1:numImgs
        %cc=conj(squeeze(img1(24,24,:,ii*2)-1*img1(24,24,:,ii*2-1)));
        delete *json *jsonhlsvd
        cc = squeeze(img1(:,ii));

        jsonmesh=struct('step_size',timeStep,'signals_i',imag((cc))'...
            ,'signals_r',real((cc))','n_singular_values',50);
        
        savejson('',jsonmesh,[prefix int2str(ii) '.json']);
        python('/vols/Scratch/asteel/forCaro/hlsvdScripts/demo_uzay.py',[prefix int2str(ii) '.json']);
        abc=loadjson([prefix int2str(ii) '.jsonhlsvd']);
        
        % tstep  = timeStep;
        % times  = 0:tstep:(512-1)*tstep;
        %  K = 1i * 2 *pi

        % for icount=1:abc.nsv_found
        % xpon   = abc.amplitudes(icount)*exp(times/abc.damping_factors(icount) + (K*( abc.frequencies(icount)*times+abc.phases(icount)/360)));
        % 
        % recon(icount,:)  = xpon' ;
        % end
        % recon(~isfinite(recon))=0;

        [~,I]=sort(abs(abc.frequencies'));

        amplitudes=abc.amplitudes(I);
        frequencies=abc.frequencies(I);
         II=find(frequencies<0.025 & frequencies > -0.025);
        damping_factors=abc.damping_factors(I);
        phases=abc.phases(I);


        tstep  = timeStep;
        times  = 0:tstep:(tps-1)*tstep;
         K = 1i * 2 *pi;
        clear recon2
        recon2 = zeros(max(II),tps);
        for icount=1:max(II)
        xpon   = amplitudes(icount)*exp(times/damping_factors(icount) +...
            (K*( frequencies(icount)*times+phases(icount)/360)));

        recon2(icount,:)  = xpon' ;

        end
        recon2(~isfinite(recon2))=0;

        %plot(real(fftshift(fft(conj(cc'-sum(recon2,1))))));
        %hold on;
        reconOut(:,ii) = conj(cc'-sum(recon2,1));
end

toc


