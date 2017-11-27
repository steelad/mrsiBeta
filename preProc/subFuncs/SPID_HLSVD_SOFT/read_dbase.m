function [stof,abs_scl_factor,filename_2fit] = read_dbase(nrp,filtertype,truncation,truncationtype)
% read database from file

pname='testdat/philips/te23/';
fname='LEEN_NORMAAL.sdat';
pnamewater='testdat/philips/te23/';
fnamewater='LEEN_WATER.sdat';
%dbasestring='PHILIPS_PRESS_TE23';
if strcmp(filtertype,'hlsvd')
dbasestring='filteredHLSVDPRO50_060107';
% if truncation
% dbasestring='filteredHLSVDPRO50_984_060107';
% end
else
dbasestring='paper051007_HLSvsFIR_PhilipsPRESS23';
end
dir_dbase_ecc = 'dbase/ecc/';
dbase_file_lip1 = 'dbase/mobile_lipids/lip1';
dbase_file_lip2 = 'dbase/mobile_lipids/lip2';

% read dbase

load([dir_dbase_ecc dbasestring])

% This loop will look for the dbase components that will be used.

nrstofs = 0;

for i = 1:8 %23
  nrstofs = nrstofs+1;
  
  % changed by JB Poullet 20060111
  if strcmp(truncationtype,'begin')&truncation 
      txttrunc = 'length(DBase)-nrp+1:length(DBase)';
  else
      txttrunc = '1:nrp';
  end
  switch dbasestring
  case 'PHILIPS_PRESS_TE23'
      signal = philips_press_te23(i,eval(txttrunc));   
  case 'PHILIPS_PRESS_TE30'
      signal = philips_press_te30(i,eval(txttrunc));
  case 'SIEMENS_STEAM_TE20'
      signal = siemens_steam_te20(i,eval(txttrunc));
  case 'SMIS_STEAM_TE10_7T'
      signal = smis_steam_te10_7t(i,eval(txttrunc));
  case 'paper051007_HLSvsFIR_PhilipsPRESS23';
      signal = DBase(i,eval(txttrunc));  
  case 'filteredHLSVDPRO50_060107';
      signal = DBase(i,eval(txttrunc));  
  case 'filteredHLSVDPRO50_984_060107';
      signal = DBase(i,eval(txttrunc));   
  end
  filename_2fit{nrstofs} = filename{i};
  stof(1:nrp,nrstofs) = signal.';
  
  % abs_scl_factor will convert results to absolute concentrations
  % variable is defined when the dbase was created.
  temp_abs_scl_factor(nrstofs,1) = abs_scl_factor(i);

end

abs_scl_factor = temp_abs_scl_factor;

% adding the mobile lipid components at 1.3ppm
load([dbase_file_lip1]) 
nrstofs = nrstofs+1;
stof(1:nrp, nrstofs)    = 0.01.*signal(1:nrp).';
filename_2fit{nrstofs}  = 'lip1';
abs_scl_factor(nrstofs,1) = 1;

% adding the mobile lipid components at 0.9ppm
load([dbase_file_lip2])
nrstofs = nrstofs+1;
stof(1:nrp, nrstofs)    = 0.01.*signal(1:nrp).';
filename_2fit{nrstofs}  = 'lip2';
abs_scl_factor(nrstofs,1) = 1;

% end read dbase