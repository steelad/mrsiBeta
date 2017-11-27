
function fidcor=phase_max_uzay(fid)
% function function fidcor=phase_max(fid);
%
% phase correction:  maximizes the real part of spectrum
% from aph0_FDmax by Pat Bolan
% A simplifid version...
% Dinesh Deelchand,  CMRR, Univ. of Minnesota

% FT of FID
spectdata = fftshift(fft(fid),1);
[np,nt]=size(spectdata);
%for kk=1:nt
 %   spectdata(:,kk)=fftshift(fft(fid(:,kk)));
%end

fidcor = complex(zeros(np,nt));
rp_out = zeros(1, nt);

for jdx = 1:nt
    spSingle = spectdata(:,jdx);
    
    % Pass #1: find the global max. Not necessary if a rp was passed in
    rps = -180:10:180;
    realsumsq = zeros(1,size(rps,2));
    for idx = 1:size(rps,2)
        spec = real(spSingle.*exp(-1i*((pi/180)*rps(idx))));
        realsumsq(idx) = sum(spec);
    end
    rpstart = rps(find(realsumsq==max(realsumsq)));
    % There can be two of these. I don't know how to choose.
    rpstart = rpstart(1);
    
    
    % Pass #2: Find the phase to 1-degree accuracy. This time use Sum of
    % squares
    clear rps realsumsq
    rps = rpstart-10:1:rpstart+10;
    realsumsq = zeros(1,size(rps,2));    
    for idx = 1:size(rps,2)
        spec= real(spSingle.*exp(-1i*((pi/180)*rps(idx))));
        realsumsq(idx) = sum(spec);
    end
    rp = rps(find(realsumsq==max(realsumsq)));
        
    % Pass #3: 0.1 degree accuracy.
    clear rps realsumsq;
    rps = rp-1:.1:rp+1;
    realsumsq = zeros(1,size(rps,2));    
    for idx = 1:size(rps,2)
        spec= real(spSingle.*exp(-1i*((pi/180)*rps(idx))));
        realsumsq(idx) = sum(spec);
    end
    rp = rps(find(realsumsq==max(realsumsq)));
    
    % Pass #4: 0.01 degree accuracy.
    clear rps realsumsq;
    rps = rp-.1:.01:rp+.1;
    realsumsq = zeros(1,size(rps,2));    
    for idx = 1:size(rps,2)
        spec= real(spSingle.*exp(-1i*((pi/180)*rps(idx))));
        realsumsq(idx) = sum(spec);
    end
    rp = rps(find(realsumsq==max(realsumsq)));
    rp_out(jdx) = rp(1);
    
    
    %Apply phase shift in Freq domain
    spectdataCorTEMP = spSingle.*exp(-1i*((pi/180)*rp_out(jdx)));
    
    % Convert back to fid.
    fidcor(:,jdx) = ifft(fftshift(spectdataCorTEMP));
end
return
