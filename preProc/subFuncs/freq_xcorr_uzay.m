function fidCor = freq_xcorr_uzay(fid,sw,sfrq1H)

% function fidCor = freq_xcorr(fid)
%
% B0 correction using cross-correlation (based on Pat Bolan freqAlignXCorr.m)
% Dinesh Deelchand, CMRR, Univ. of Minnesota
% 31 March 2011

%global sw

% define shift threshold in Hz
maxshiftHz = 100;

%%%%%%%%%%%%%%%%%%%%%%%%
% apply LB and ZF
dw = 1/sw; t = (0:dw:dw*(length(fid)-1))';
LB = 10;
GF = 1000;
sifactor = 10;
[np,nt,nbCoils]=size(fid);
fidzf = complex(zeros(np*(sifactor+1),nt,nbCoils));
fidCor = complex(zeros(np,nt,nbCoils));
for ical=1:nbCoils
    for jcal=1:nt
        fidzf(:,jcal,ical) = [fid(:,jcal,ical).*exp(-t*pi*LB-t.^2/(GF^2)); zeros(np*sifactor,1)];
    end
end

%h = waitbar(0,'Running cross-correlation algorithm, please wait...');

%%%%%%%%%%%%%%%%%%%%%%%%
% FFT and mean
spectfft = fftshift(fft(fidzf,[],1),1);
spect_ref = abs(mean(mean(spectfft,2),3));

%%%%%%%%%%%%%%%%%%%%%%%%
% Cross-correlation algorithm
fmax = (sw)/2;
f = fmax:-2*fmax/(length(fidzf)-1):-fmax;
newSR = sw./length(spect_ref);
maxlags = round(maxshiftHz./newSR);
FreqRef = 0;
count = 0;
freqC = zeros(nt,nbCoils);
diffFreq = zeros(nt,nbCoils);
for iy = 1:nbCoils
   % fprintf('     %2g ',iy)
    for ix = 1:nt
        clear spect_use
        spect_use = abs(spectfft(:,ix,iy));
        c = xcorr(spect_ref, spect_use, maxlags);
        
        % Determine frequency shift based on position of max signal
        [~,inx] = max(c);
        freqC(iy,ix) = f(inx);
        
        % Frequency shift (in Hz) wrt to first data point
        if (ix==1) && (iy==1)
            FreqRef = freqC(1,1);
        end
        diffFreq(iy,ix) = freqC(iy,ix) - FreqRef;
        
        % Correct FIDs with determined frequency shift
        fidCor(:,ix,iy) = fid(:,ix,iy) .*exp(1i*2*pi*-diffFreq(iy,ix)*t);
        
        % waitbar
        count = count + 1;
       % waitbar(count/(nt*nbCoils));
    end
end
%close(h);
%fprintf('\n')
plotflag=1;
H1offset=4.67;
if plotflag==1
    % apply LB and ZF
    dw = 1/sw; t = (0:dw:dw*(length(fid)-1))';
    LB = 5; GF = 0.15; sifactor = 3;
    [np,nt,nbCoils] = size(fid);
    fidzf = complex(zeros(np*(sifactor+1),nt,nbCoils));
    fidzfCor = fidzf;
    fmax = (sw)/2;
    f = fmax:-2*fmax/(length(fidzf)-1):-fmax;
    scale_ppm = f/(sfrq1H)+H1offset;

    for ical = 1:nbCoils
        for jcal = 1:nt
            fidzf(:,jcal,ical) = [fid(:,jcal,ical).*exp(-t*pi*LB-t.^2/(GF^2)); zeros(np*sifactor,1)];
            fidzfCor(:,jcal,ical) = [fidCor(:,jcal,ical).*exp(-t*pi*LB-t.^2/(GF^2)); zeros(np*sifactor,1)];
        end
    end
   % figure;
    
        spectfftold = real(fftshift(fft(fidzf,[],1),1));
   
    %subplot(211), plot(scale_ppm,real(spectfftold(:,:))), title('Before Frequency correction')
    %set(gca,'xdir','reverse'); grid on;
    %curaxis=axis; axis([0 5 curaxis(3) curaxis(4)]); newaxis=axis;
    

        spectfftnew = real(fftshift(fft(fidzfCor,[],1),1));
        
  
    %subplot(212), plot(scale_ppm,spectfftnew(:,:)), title('After correction');
    %set(gca,'xdir','reverse'); grid on;
    %axis(newaxis);
end

return