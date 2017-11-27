% baseline construction

close all
clear all

%Paper by Seeger Magn Reson Med 49:19-28 (2003)
% Gaussian peaks in 1.51,2.05,2.1,2.24,2.81,3.0

dim = 1024; %number of points
nbbas = 1; %number of baselines 
t = [0:dim-1];
frequency = 63.898; %spectrometer frequency in MHz
ngc = 1; %number of gauss components 


%meang = [1.51 2.05 2.1 2.24 2.81 3.0];
% Position: mean value - in ppm
meanmac = [1.51];% 2.1 3.0];
meanlip = [2.05];% 2.24 2.81];

% meang = frequency*(meang-4.7); %in Hz
%sigmag = [39.5 13.1 27.5 16.5 16.1 14.0]/2;

% Linewidths (LW) - in Hz (width at midheight)
LWmac = [25];% 27.5 14.0]; %macromolecules
LWlip = [13.1];% 16.5 16.1]; %lipids (except lipids at 1.3 and 0.89 ppm)

% standard deviation of the gaussian - in Hz (sigma gauss components) 
% density function of a normal distribution centered in 0:
% f(x)=((sigma*sqrt(2*pi))^-1)*exp(x^2/(2*sigma^2))
% LW = 2*x = 2*sqrt(-2*sigma^2*ln(0.5))
% sigma = sqrt(-x^2/(2*ln(0.5)))
sigmac = sqrt(-(LWmac/2).^2/(2*log(0.5)));
siglip = sqrt(-(LWlip/2).^2/(2*log(0.5)));

%  sigmag = [39.5 13.1 27.5 16.5 16.1 14.0]; % in Hz 

% In ppm
%   sigmac = sigmac/frequency; %in ppm
%   siglip = siglip/frequency; %in ppm

 % ampl_gauss_components   = [0.92 0.15 0.76 0.16 0.10 0.11];
  ampl_gauss_components_mac   = [0.92];% 0.76 0.11];
  ampl_gauss_components_lip   = [0.15];% 0.16 0.10];

[xpts,xpts2] = absi(dim-1,1); % [in kHz,in ppm]

baselinemac = zeros(1,dim); 
baselinelip = zeros(1,dim); 

% Macromolecules
for nos=1:nbbas
  agcmac = ampl_gauss_components_mac;%+ampl_gauss_components/10/1.96.*randn(1,length(ampl_gauss_components));
  for i = 1:ngc
   %baselinemac = baselinemac + agcmac(i)*normpdf(xpts2,meanmac(i),sigmac(i)); %frequency domain (real part)
   baselinemac = baselinemac + agcmac(i)*exp(-(LWmac(i)*sqrt(2))^2*(t/1000).^2+sqrt(-1)*(2*pi*(meanmac(i)-4.7)*frequency*t/1000));
  %baselinemac = baselinemac + agcmac(i)*1/sqrt(sigmac(i)^2*pi)*exp(-xpts2/sigmac(i)^2).*exp(sqrt(-1)*2*pi*meanmac(i)*xpts2);
  end
end %for nos=1:nbbas

% basemac = zeros(1,dim)+ agcmac(i)*normpdf(xpts2,meanmac(i),sigmac(i));
% basemac= ifft(ifftshift(basemac));
% figure
% disp_sig2(basemac,0,1,0,1024,'',0,1,'nmr',1,0);

% for nos=1:nbbas
%   agclip = ampl_gauss_components_lip;%+ampl_gauss_components/10/1.96.*randn(1,length(ampl_gauss_components));
%   for i = 1:ngc
%    baselinelip = baselinelip + agclip(i)*normpdf(xpts2,meanlip(i),siglip(i)); %frequency domain (real part)
%   end
% end %for nos=1:nbbas
%baselinelip = zeros(1,length(baselinemac));

baseline = baselinemac; %+ baselinelip;

%baseline = baseline +baseline*sqrt(-1);

%DISPLAYING
%baseline= ifft(ifftshift(baseline));
% % length(baseline)/2
% % baseline = [baseline(1:length(baseline)/2) zeros(1,length(baseline)/2)];
% baselinemac= ifft(ifftshift(baselinemac));
% baselinelip= ifft(ifftshift(baselinelip));

figure
plot(real(baseline));

% figure
% plot3(real(baseline),imag(baseline),[0:dim-1]);

freq = frequency*1000;

figure
disp_sig2(baseline,0,1,0,1024,'',0,0,freq,0,0);

figure
disp_sig2(baseline,0,1,0,1024,'',0,2,freq,0,0);

% figure
% disp_sig2(baselinemac,0,1,0,1024,'',0,1,'nmr',0,0);

% length(baseline)/2
% baseline = [baseline(1:length(baseline)/2) zeros(1,length(baseline)/2)];
% figure
% plot(real(baseline));
% 
% figure
% disp_sig2(baseline,0,1,0,1024,'',0,1,'nmr',0,0);
% hold on
% disp_sig2(2*baseline,0,1,0,1024,'',0,1,'nmr',0,0,'r');


% figure
% disp_sig2(baselinelip,0,1,0,1024,'',0,1,'nmr',1,0);