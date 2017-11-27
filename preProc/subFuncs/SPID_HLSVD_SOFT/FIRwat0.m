%FIR filtering the "Water" components

function h=FIRwat0(PB,signal,order,step,ndp,TBW,Dfact,dis)
%signal = signal to be filtered
%step = time step
%ndp = number od data points
%PB = passband in kHz (frequency )
%order = filter order
%disp = 1 if you want to display results
%TBW = transition band width (in kHz), typical value is 0.02
%Dfact = delay factor
% h = filter coefficients

if size(signal,1)>size(signal,2)
    signal=signal.'; %transposed (non conjugated)
end

%% Definition frequency and damping grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<8
    dis=0;
end
if nargin<7
    Dfact=2.2/3;
end
if nargin<6
    TBW=0.02;
end

bestminR = [];
besth=[];
bestfgs=[];
bestfgp=[];
bestH1=[];
besth1=[];
bestdg=[];
bestpM2=[];
bestN=[];



N = order+1;                % numerator order discrete-time filter
fpu0 =  (PB(2)-PB(1))/2;
fpl0 = -fpu0;
fsu0 = fpu0+TBW;
fsl0 = -fsu0;

fsl02 = PB(1)-TBW;
% fpl02 = PB(1);
% fpu02 = PB(2);
fsu02 = PB(2)+TBW;
D = round(N*Dfact); %delay
j=sqrt(-1);

noise = std(signal(end-20:end))/sqrt(2); %std of the noise in the time domain
noise =  noise /sqrt(ndp); %std of the noise in the frequency domain
if noise <1e-6
    noise=1e-6;
end

fmin = -0.5;
fmax = 0.5;
kHz = -0.5:1/length(signal):0.5-1/length(signal);

%% peak picking
sf = fftshift(fft(signal));

[m,indx]=find(fsl02>=kHz | fsu02<=kHz);
indx(indx>length(sf))=[];
%         figure
%         plot(abs(sf(indx)))
[md,id]=max(abs(sf(indx)));

f1 = [];

dmin = 0;
dmax = 0.03;

minRold = abs(md);
bestminR =minRold;
besth=[];
iiv=1;
c=1; %counting the iterations
while iiv>0
    fw=kHz(indx(id+1));
    if (fw>fsl02) && (fw<fsu02)
        warning('The water resonances should be in the stop band, not in the transition band. The design of the filter has been interrupted.')
        return
    end
    if fw > fsu02
        %p=1; %frequency larger that the frequencies of the passband
        fw2 = fsl02-(fw-fsu02);
    else
        %p=0; %frequency smaller that the frequencies of the passband
        fw2 = fsu02+(fsl02-fw);
    end
    [fwT]=sort([fw fw2]); %fw = lower frequency water resonance, fw2  = higher frequency water resonance

    fw(2)=(fwT(2)-fwT(1))/2;
    fw(1)=-fw(2);
%     while sum(f1==fw(2))>0 % to avoid having several zeroes on the same z-plane location.
%         fw(2) = fw(2)-0.001;
%         fw(1) = -fw(2);
%     end
    f1 = [f1 fw];
    z1 = exp(j*2*pi*f1);
    q1 = 1./z1;

    if dis
        figure
        plot(z1,'ro')
        set(gca,'xlim',[-1,1])
        hold on
        o = -pi:pi/30:pi;
        x = cos(o);
        y = sin(o);
        plot(x,y)
        title('Roots chosen for H1')
    end
    %%
    %%%%%%%%%%
    % Choosing Grid
    %%%%%%%%%%
    fgp = fpl0:0.001:fpu0; % only passband
    f1 = sort(f1);
    fgs = [-0.5:0.001:fsl0 fsu0:0.001:0.5];
    ind = [];
    for i=1:length(fgs)
        if min(abs(fgs(i)-f1))<0.001
            ind = [ind i];
        end
    end
    fgs(ind) = [];
    fg = [fgp fgs];

    if c==1;
        dg = 0;
    end

    %z2 = exp(sqrt(-1)*2*pi*fg');% * exp(-dg);
    z2 = exp(sqrt(-1)*2*pi*fg')*exp(-dg);
    zt2 = log(z2);
    zM2 = exp(zt2(:)*(0:N-1-length(q1)));
    qM2 = 1./zM2;
    q2 = 1./z2;

    %%
    %%%%%%%%%%
    % H1 computation
    %%%%%%%%%%
    H1 = ones(size(z2,1)*size(z2,2),1);
    for i=1:length(z2)
        for ii=1:length(f1)
            H1(i) = H1(i)*(1-q1(ii)*q2(i));
        end
    end

    %%
    %%%%%%%%%%
    % H2 computation
    %%%%%%%%%%
    TT = repmat([(exp(-sqrt(-1)*D*2*pi*fgp'));zeros(length(fgs),1)],1,length(dg));
    T = TT(:)./H1;
    W = abs(H1).^2;

    pM2 = qM2;
    %h2 = pM2\T;
    %cond(pM2)
    h2 = (pM2'*diag(W)*pM2) \ (pM2'*diag(W)*T);

    %%
    %%%%%%%%%%
    % H=H1*H2
    %%%%%%%%%%
    h1 = poly(leja(q1));
    h = conv(h1,h2);
    %h=h2;
    h=minphJB(h);
    %h=h';
    wc = PB(1)+(PB(2)-PB(1))/2;
    h = h.*exp(j*2*pi*(wc)*[0:N-1]);

    %%
    %%%%%%%%%%
    % plotting results: |H| wrt f and d
    %%%%%%%%%%
    if dis
        figure
        fg5 = (fmin:(fmax-fmin)/100:fmax)';
        dg5 = dmin:(dmax-dmin)/100:dmax;
        z5 = exp(sqrt(-1)*2*pi*fg5) * exp(-dg5);
        zt5 = log(z5);
        zM5 = exp(zt5(:)*(0:N-1));
        q5 = 1./z5(:);
        qM5 = 1./zM5;
        H5 = qM5(:,1:length(h))*h(:);
        fgl5 = length(fg5);
        dgl5 = length(dg5);
        mesh(dg5,fg5,abs(reshape(H5,fgl5,dgl5)))
        xlabel('damping')
        ylabel('normalized frequency')
        zlabel('freq-damp response')

        figure
        plot(fg5,db(abs(H5(1:length(fg5)))))
        xlabel('normalized frequency')
        ylabel('frequency response')

        % plotting filter caracteristics
        %         fvtool(h1)
        %         fvtool(h2)
        %         fvtool(h)
    end
    h = fliplr(h);
    h = conj(h);
    sf=filter(h,1,signal);
    sf=[sf(length(h)-1:end) zeros(1,length(h)-2)];
    sfF = fftshift(fft(sf));
    dg = findDamping(sf,step,ndp);
    if dg>0.015
        dg=0.015;
    end

    if dis
        fq = [-0.5:1/length(signal):0.5-1/length(signal)];
        figure
        plot(fq,abs(sfF),'b')
        hold on
        plot(fq,abs(fftshift(fft(signal))),'k')
    end

    %% Checking the errors
    if id+1>length(indx)
        [md,id]=max(abs(sfF(indx(1:end-1))));
        %disp(['The maximum residual in the stop band is in fact' num2str(md)]);
    else
        [md,id]=max(abs(sfF(indx)));
        %disp(['The maximum residual in the stop band is ' num2str(md)]);
    end
    minR=abs(md)/sqrt(length(signal));
%     disp(['The filter order is ',num2str(N)])
%     disp(['The zeroes are  ',num2str(f1)])
%     disp(['The 2.5*noise is ',num2str(2.5*noise)])
%     disp(['minR ',num2str(minR)])
%     disp(['minRold ',num2str(minRold)])

    if minR<bestminR
        %checking results to avoid recording bad filter coefficients
        fg5 = (fmin:(fmax-fmin)/1000:fmax)';
        dg5 = dmin:(dmax-dmin)/5:dmax;
        z5 = exp(sqrt(-1)*2*pi*fg5) * exp(-dg5);
        zt5 = log(z5);
        zM5 = exp(zt5(:)*(0:N-1));
        q5 = 1./z5(:);
        qM5 = 1./zM5;
        H5 = qM5(:,1:length(h))*h(:);
        [mF,indxF]=find(fsl02<fg5' & fsu02>fg5');
        indL=round(length(indxF)/2);
        diff=abs(sum(abs(H5(indxF(round(indL/2):indL))))-sum(abs(H5(indxF(indL+1:end-round(indL/2)+1)))));
        %diff=abs(sum(abs(H5(indxF(1:indL))))-sum(abs(H5(indxF(indL+1:end)))))
        if diff<10
            bestminR = minR;
            besth=h;
            bestfgs=fgs;
            bestfgp=fgp;
            bestH1=H1;
            besth1=h1;
            bestdg=dg;
            bestpM2=pM2;
            bestN=N;
        end
    end

    if minR<2.5*noise
        z2 = exp(sqrt(-1)*2*pi*fg')*exp(-dg);
        zt2 = log(z2);
        zM2 = exp(zt2(:)*(0:N-1-length(q1)));
        qM2 = 1./zM2;
        q2 = 1./z2;

        %%
        %%%%%%%%%%
        % H1 computation
        %%%%%%%%%%
        H1 = ones(size(z2,1)*size(z2,2),1);
        for i=1:length(z2)
            for ii=1:length(f1)
                H1(i) = H1(i)*(1-q1(ii)*q2(i));
            end
        end

        %%
        %%%%%%%%%%
        % H2 computation
        %%%%%%%%%%
        TT = repmat([(exp(-sqrt(-1)*D*2*pi*fgp'));zeros(length(fgs),1)],1,length(dg));
        T = TT(:)./H1;
        W = abs(H1).^2;

        pM2 = qM2;
        %h2 = pM2\T;
        %cond(pM2)
        h2 = (pM2'*diag(W)*pM2) \ (pM2'*diag(W)*T);

        %%
        %%%%%%%%%%
        % H=H1*H2
        %%%%%%%%%%
        %h1 = poly(leja(q1));
        h = conv(h1,h2);
        %h=h2;
        h=minphJB(h);
        %h=h';
        wc = PB(1)+(PB(2)-PB(1))/2;
        h = h.*exp(j*2*pi*(wc)*[0:N-1]);
        h = fliplr(h);
        h = conj(h);

        %checking results
        fg5 = (fmin:(fmax-fmin)/1000:fmax)';
        dg5 = dmin:(dmax-dmin)/5:dmax;
        z5 = exp(sqrt(-1)*2*pi*fg5) * exp(-dg5);
        zt5 = log(z5);
        zM5 = exp(zt5(:)*(0:N-1));
        q5 = 1./z5(:);
        qM5 = 1./zM5;
        H5 = qM5(:,1:length(h))*h(:);
        [mF,indxF]=find(fsl02<fg5' & fsu02>fg5');
%         figure
%         plot(db(abs(H5(indxF))))

        indL=round(length(indxF)/2);
        diff=abs(sum(abs(H5(indxF(round(indL/2):indL))))-sum(abs(H5(indxF(indL+1:end-round(indL/2)+1)))));
        bestdiff=diff;
        %optimizing the delay
        bt=1;
        DfactT = Dfact;
        prev = [];
        prevnew=[];
        bestminR = [];
        if isempty(besth)
            besth=h;
            bestfgs=fgs;
            bestfgp =fgp;
            bestH1 = H1;
            besth1= h1;
            bestdg = dg;
            bestpM2 = pM2;
            bestN=N;
        end
        
        DS = 5;
        initstep = min(1-DfactT,DfactT-0);
        cc=0;
        while bt==1
            DfactTS = initstep/DS; %Dfact step
            %init Dfacts
            DfactTP = [DfactT-DfactTS,DfactT+DfactTS];
            adbt=0;
            for it=1:2
                %D = round(bestN*DfactTP(it)); %delay
                D = bestN*DfactTP(it); %delay
                %%%%%%%%%%
                % H2 computation
                %%%%%%%%%%
                TT = repmat([(exp(-sqrt(-1)*D*2*pi*bestfgp'));zeros(length(bestfgs),1)],1,length(bestdg));
                T = TT(:)./bestH1;
                W = abs(bestH1).^2;
                h2 = (bestpM2'*diag(W)*bestpM2) \ (bestpM2'*diag(W)*T);
                %%%%%%%%%%
                % H=H1*H2
                %%%%%%%%%%
                h = conv(besth1,h2);
                h=minphJB(h);
                wc = PB(1)+(PB(2)-PB(1))/2;
                h = h.*exp(j*2*pi*(wc)*[0:bestN-1]);
                h = fliplr(h);
                h = conj(h);
                H5 = qM5(:,1:length(h))*h(:);
                indL=round(length(indxF)/2);
                diff=abs(sum(abs(H5(indxF(round(indL/2):indL))))-sum(abs(H5(indxF(indL+1:end-round(indL/2)+1)))));
                if bestdiff>diff
                    bestdiff=diff;
                    prevnew = it;
                    DfactT = DfactTP(it);
                    adbt=adbt+1;
                    besth=h;
                end
            end
            if prev ~= prevnew
                DS = DS*2;
            else
                DS = DS*1.1;
            end
            if adbt<1
                     cc=cc+1;
                    DS = DS*2;
                    if cc==5
                        bt=0;
                    end
             end
        end
        h=besth;
        diff=abs(sum(abs(H5(indxF(round(indL/2):indL))))-sum(abs(H5(indxF(indL+1:end-round(indL/2)+1)))));
        H5 = qM5(:,1:length(h))*h(:);
%         figure
%         plot(db(abs(H5(indxF))))
        
        if dis
        figure
        fgl5 = length(fg5);
        dgl5 = length(dg5);
        mesh(dg5,fg5,abs(reshape(H5,fgl5,dgl5)))
        xlabel('damping')
        ylabel('normalized frequency')
        zlabel('freq-damp response')
        end
        
        iiv=0;
    end
    if minRold<minR
        N = N+10; % changing the filter order
        if N>110
            warning('The attenuation conditions have not been met. The chosen filter may not be the optimal filter.')
            %we keep the best h, i.e. the one providing the smallest minR
         if isempty(besth)
            besth=h;
            bestfgs=fgs;
            bestfgp =fgp;
            bestH1 = H1;
            besth1= h1;
            bestdg = dg;
            bestpM2 = pM2;
            bestN=N;
         else
            h=besth;
         end
%             if dis
%                 %checking results
%                 fg5 = (fmin:(fmax-fmin)/100:fmax)';
%                 dg5 = dmin:(dmax-dmin)/100:dmax;
%                 z5 = exp(sqrt(-1)*2*pi*fg5) * exp(-dg5);
%                 zt5 = log(z5);
%                 zM5 = exp(zt5(:)*(0:N-1));
%                 q5 = 1./z5(:);
%                 qM5 = 1./zM5;
%                 H5 = qM5(:,1:length(h))*h(:);
%                 fgl5 = length(fg5);
%                 dgl5 = length(dg5);
%                 figure
%                 mesh(dg5,fg5,abs(reshape(H5,fgl5,dgl5)))
%                 xlabel('damping')
%                 ylabel('normalized frequency')
%                 zlabel('freq-damp response')
%             end
                        %checking results
            fg5 = (fmin:(fmax-fmin)/1000:fmax)';
            dg5 = dmin:(dmax-dmin)/5:dmax;
            z5 = exp(sqrt(-1)*2*pi*fg5) * exp(-dg5);
            zt5 = log(z5);
            zM5 = exp(zt5(:)*(0:N-1));
            q5 = 1./z5(:);
            qM5 = 1./zM5;
            H5 = qM5(:,1:length(h))*h(:);
            [mF,indxF]=find(fsl02<fg5' & fsu02>fg5');
%             figure
%             plot(db(abs(H5(indxF))))

            indL=round(length(indxF)/2);
            diff=abs(sum(abs(H5(indxF(round(indL/2):indL))))-sum(abs(H5(indxF(indL+1:end-round(indL/2)+1)))));
            %diff=abs(sum(abs(H5(indxF(1:indL))))-sum(abs(H5(indxF(indL+1:end)))));
            bestdiff=diff;
            %optimizing the delay
            bt=1;
            DfactT = Dfact;
            prev = [];
            prevnew=[];
            DS = 5;
            cc=0;
            initstep = min(1-DfactT,DfactT-0);
            while bt==1
                DfactTS = initstep/DS; %Dfact step
                %init Dfacts
                DfactTP = [DfactT-DfactTS,DfactT+DfactTS];
                adbt=0;
                for it=1:2
                    %D = round(bestN*DfactTP(it)); %delay
                    D = bestN*DfactTP(it); %delay
                    %%%%%%%%%%
                    % H2 computation
                    %%%%%%%%%%
                    TT = repmat([(exp(-sqrt(-1)*D*2*pi*bestfgp'));zeros(length(bestfgs),1)],1,length(bestdg));
                    T = TT(:)./bestH1;
                    W = abs(bestH1).^2;
                    h2 = (bestpM2'*diag(W)*bestpM2) \ (bestpM2'*diag(W)*T);
                    %%%%%%%%%%
                    % H=H1*H2
                    %%%%%%%%%%
                    h = conv(besth1,h2);
                    h=minphJB(h);
                    wc = PB(1)+(PB(2)-PB(1))/2;
                    h = h.*exp(j*2*pi*(wc)*[0:bestN-1]);
                    h = fliplr(h);
                    h = conj(h);
                    H5 = qM5(:,1:length(h))*h(:);
                    indL=round(length(indxF)/2);
                    diff=abs(sum(abs(H5(indxF(round(indL/2):indL))))-sum(abs(H5(indxF(indL+1:end-round(indL/2)+1)))));
                    %diff=abs(sum(abs(H5(indxF(1:indL))))-sum(abs(H5(indxF(indL+1:end)))));
                    if bestdiff>diff
                        bestdiff=diff;
                        prevnew = it;
                        DfactT = DfactTP(it);
                        adbt=adbt+1;
                        besth=h;
                    end
                end
                if prev ~= prevnew
                    DS = DS*2;
                else
                    DS = DS*1.1;
                end
                if adbt<1
                    cc=cc+1;
                    DS = DS*2;
                    if cc==5
                        bt=0;
                    end
                end
            end
            h=besth;
%             H5 = qM5(:,1:length(h))*h(:);
%             figure
%             fgl5 = length(fg5);
%             dgl5 = length(dg5);
%             mesh(dg5,fg5,abs(reshape(H5,fgl5,dgl5)))
%             xlabel('damping')
%             ylabel('normalized frequency')
%             zlabel('freq-damp response')
            %plot(db(abs(H5(indxF))))
            return
        end
        sf = fftshift(fft(signal));
        %[m,indx]=find(fsl02>=kHz | fsu02<=kHz);
        %indx(indx>length(sf))=[];
        [md,id]=max(abs(sf(indx)));
        f1 = [];
        minRold = abs(md);
        %D = round(N*Dfact); %delay
        D = N*Dfact; %delay
    else
        minRold = minR;
    end
    c=c+1;
    %    disp('----------------------------------------')
end

