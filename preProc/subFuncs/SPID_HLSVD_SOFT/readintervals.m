function [compNames,freqmatrix] = readintervals(ir)

if ir ==1 % short-echo time MR spectra at 1.5 T
 compNames = {'Ala' ;'Gaba'; 'NAA'; 'Gln' ;...
        'Gaba2' ;'Cr' ;'Cho' ;'PCh ';'Tau ';'Myo';'Glu';'Cr1';'Lac'};
    peakwidth = 0.04;
    % center frequency
    Ala = 1.48;
    Gaba = 1.9;
    NAA = 2.03;
    Gln = 2.11;
    Gaba2 = 2.27;
    Cr = 3.03;
    Cho = 3.19;
    PCh = 3.215;
    GPCh = 3.239;
    Tau = 3.42;
    Myo = 3.54;
    Glu = 3.76;
    Cr1 = 3.93;
    Lac = 4.12;

    peaks(1,:) = [1.46 1.50]; % Ala
    peaks(2,:) = [1.875 1.925]; % Gaba
    peaks(3,:) = [1.99  2.025]; % NAA    
    peaks(4,:) = [2.09 2.17]; % Glu (+ Gln)
    peaks(5,:) = [2.245 2.295]; % Gaba2
    peaks(6,:) = [3.01 3.05]; % Cr
    peaks(7,:) = [3.19 3.205]; % Cho
    peaks(8,:) = [3.205 3.245]; % PCh
    peaks(9,:) =[3.40 3.45]; % Tau
    peaks(10,:) =[3.50 3.58]; % Myo
    peaks(11,:) =[3.71 3.82]; % Glu3
    peaks(12,:) =[3.91 3.945]; % Cr2
    peaks(13,:) =[4.08 4.14]; % Lac2
    
    freqmatrix = peaks;
end


if ir ==2 % long echo time in vivo MR spectra 
 compNames = {'Ala' ;'Gaba'; 'NAA'; 'Gln' ;...
        'Gaba2' ;'Cr' ;'Cho' ;'PCh ''Tau ';'Myo';'Glu';'Cr1';'Lac'};
    peakwidth = 0.04;
    % center frequency
    Ala = 1.48;
    Gaba = 1.9;
    NAA = 2.03;
    Gln = 2.11;
    Gaba2 = 2.27;
    Cr = 3.03;
    Cho = 3.19;
    PCh = 3.215;
    GPCh = 3.239;
    Tau = 3.42;
    Myo = 3.54;
    Glu = 3.76;
    Cr1 = 3.93;
    Lac = 4.12;

    peaks(1,:) = [1.46 1.50]; % Ala
    peaks(2,:) = [1.875 1.925]; % Gaba
    peaks(3,:) = [1.99  2.025]; % NAA    
    peaks(4,:) = [2.09 2.17]; % Glu (+ Gln)
    peaks(5,:) = [2.245 2.295]; % Gaba2
    peaks(6,:) = [3.01 3.05]; % Cr
    peaks(7,:) = [3.19 3.205]; % Cho
    peaks(8,:) = [3.205 3.245]; % PCh
    peaks(9,:) =[3.40 3.45]; % Tau
    peaks(10,:) =[3.50 3.58]; % Myo
    peaks(11,:) =[3.71 3.82]; % Glu3
    peaks(12,:) =[3.91 3.945]; % Cr2
    peaks(13,:) =[4.08 4.14]; % Lac2
         
    freqmatrix = peaks;
end


if ir ==3 % ex vivo HRMAS MR spectra at 500 MHz
  compNames = {'Lac';'Ala' ;'Gaba'; 'NAA'; 'Glu' ;...   %Ukn
         'Glu2' ;'Cr' ;'Cho' ;'PCh ';'GPCh';'Tau ';'Pro';...
        'Tau2';'Myo';'Myo2';'Glu3';'Cr2';'Myo3';'Lac2'};
    peaks(1,:) = [1.30 1.34]; % Lac 
    peaks(2,:) = [1.45 1.49]; % Ala
    peaks(3,:) = [1.84 1.94]; % Gaba
    peaks(4,:) = [1.99  2.025]; % NAA
    peaks(5,:) = [2.09 2.17]; % Glu (+ Gln)
%    peaks(6,:) = [2.21 2.23]  % Ukn
    peaks(6,:) = [2.39 2.50]; % Glu2
    peaks(7,:) = [3.01 3.03]; % Cr
    peaks(8,:) = [3.19 3.205]; % Cho
    peaks(9,:) = [3.205 3.23]; % PCh
    peaks(10,:) = [3.235 3.245]; %GPCh
    peaks(11,:) =[3.255 3.275]; % Tau
    peaks(12,:) =[3.33 3.36]; % Pro
    peaks(13,:) =[3.39 3.435]; % Tau2 
    peaks(14,:) =[3.50 3.58]; % Myo
    peaks(15,:) =[3.58 3.70]; % Myo2
    peaks(16,:) =[3.71 3.82]; % Glu3
    peaks(17,:) =[3.91 3.945]; % Cr2
    peaks(18,:) =[4.03 4.075]; % Myo3
    peaks(19,:) =[4.08 4.14]; % Lac2
    
    freqmatrix = peaks;
end

if ir == 4 % ex vivo HRMAS MR spectra at 600 MHz
      compNames = {'Lac';'Ala' ;'Gaba'; 'NAA'; 'Glu' ;...%Ukn
         'Glu2' ;'Cr' ;'Cho' ;'PCh ';'GPCh';'Tau ';'Pro';...
        'Tau2';'Myo';'Myo2';'Glu3';'Cr2';'Myo3';'Lac2'};
    peaks(1,:) = [1.30 1.34]; % Lac 
    peaks(2,:) = [1.45 1.49]; % Ala
    peaks(3,:) = [1.84 1.94]; % Gaba
    peaks(4,:) = [1.99  2.025]; % NAA
    peaks(5,:) = [2.09 2.17]; % Glu (+ Gln)
%    peaks(6,:) = [2.21 2.23]  % Ukn
    peaks(6,:) = [2.39 2.50]; % Glu2
    peaks(7,:) = [3.01 3.03]; % Cr
    peaks(8,:) = [3.19 3.205]; % Cho
    peaks(9,:) = [3.205 3.23]; % PCh
    peaks(10,:) = [3.235 3.245]; %GPCh
    peaks(11,:) =[3.255 3.275]; % Tau
    peaks(12,:) =[3.33 3.36]; % Pro
    peaks(13,:) =[3.39 3.435]; % Tau2 
    peaks(14,:) =[3.50 3.58]; % Myo
    peaks(15,:) =[3.58 3.70]; % Myo2
    peaks(16,:) =[3.71 3.82]; % Glu3
    peaks(17,:) =[3.91 3.945]; % Cr2
    peaks(18,:) =[4.03 4.075]; % Myo3
    peaks(19,:) =[4.08 4.14]; % Lac2
    
    freqmatrix = peaks;
end    
if ir>5
    error('No set of intervals corresponds to your selection. You have to choose between 1, 2, 3 or 4.')
    return
end
