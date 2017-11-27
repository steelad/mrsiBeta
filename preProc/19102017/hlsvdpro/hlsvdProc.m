function reconOut = hlsvdProc(img1,prefix)
%%

disp('++ HLSVD being applied to all images in set... This will take a while.')
%%
disp(size(img1))
numImgs = size(img1,4);
reconOut = zeros(size(img1));
ii = 1;
tic
for vX = 1:size(img1,1)
    disp(['++ Starting vX : ' num2str(vX)]);
    for vY = 1:size(img1,2)
        disp(['++ Starting vY : ' num2str(vY)]);
        for ii=1:numImgs
        %cc=conj(squeeze(img1(24,24,:,ii*2)-1*img1(24,24,:,ii*2-1)));
        delete *json *jsonhlsvd
        cc = squeeze(img1(vX,vY,:,ii));

        jsonmesh=struct('step_size',0.8,'signals_i',imag((cc))'...
            ,'signals_r',real((cc))','n_singular_values',50);

        savejson('',jsonmesh,[prefix int2str(ii) '.json']);
        python('./hlsvdpro_copy/demo_uzay.py',[prefix int2str(ii) '.json']);
        abc=loadjson([prefix int2str(ii) '.jsonhlsvd']);
        
        % tstep  = 0.8;
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


        tstep  = 0.8;
        times  = 0:tstep:(1024-1)*tstep;
         K = 1i * 2 *pi;
        clear recon2
        recon2 = zeros(max(II),1024);
        for icount=1:max(II)
        xpon   = amplitudes(icount)*exp(times/damping_factors(icount) +...
            (K*( frequencies(icount)*times+phases(icount)/360)));

        recon2(icount,:)  = xpon' ;

        end
        recon2(~isfinite(recon2))=0;

        %plot(real(fftshift(fft(conj(cc'-sum(recon2,1))))));
        %hold on;
        reconOut(vX,vY,:,ii) = conj(cc'-sum(recon2,1));
        end
    end
end
toc


