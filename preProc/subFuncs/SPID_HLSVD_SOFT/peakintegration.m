% peak integration 
function [freqmatrix, results] = peakintegration(signal, freqmatrix,axis)
spacing = axis(2)-axis(1);
%the calculation needs to be done for all the measurements in signal
for k=1:1:size(signal,1)

    for m=1:1:size(freqmatrix,1)
        %indices of the start and stop of the intergrationinterval
        freqmatrix(m,3) = find(axis<=(freqmatrix(m,1)+(spacing/2)) & axis>=(freqmatrix(m,1)-(spacing/2)));
        freqmatrix(m,4) = find(axis<=(freqmatrix(m,2)+(spacing/2)) & axis>=(freqmatrix(m,2)-(spacing/2)));

        %real start and stop of the integrationinterval
        freqmatrix(m,5) = axis(freqmatrix(m,3));
        freqmatrix(m,6) = axis(freqmatrix(m,4));

        %integrate
        %freqmatrix(m,8)=trapz(spectra(k,freqmatrix(m,4):spacing:freqmatrix(m,5)))*spacing;
        freqmatrix(m,7)=trapz(signal(k,freqmatrix(m,3):freqmatrix(m,4)))*spacing;
    end
    %store the integrationvalues 
    temp = freqmatrix(:,7);
    results(k,:) = temp';
    freqmatrix(:,7) = [];
end
