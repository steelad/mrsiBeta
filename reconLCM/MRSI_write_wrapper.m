function [out] = MRSI_write_wrapper(hdr,metab_map,vSHIFT,outname)
%%

[~,other]    = select_tomographic_images(hdr);
[spect,~]      = select_spectroscopy_images(other);
fspe = convert_spectroscopy(spect,'flat','nii',metab_map,vSHIFT,outname);


out.files = [fspe(:)];
if isempty(out.files)
    out.files = {''};
end;
return;
%%
function fnames = convert_spectroscopy(hdr,root_dir,format,metab_map,vSHIFT,outname)
fnames = cell(length(hdr),1);
for i=1:length(hdr),
    fnames{i} = write_spectroscopy_volume(hdr(i),root_dir,format,metab_map,vSHIFT,outname);
end;
return;

function fname = write_spectroscopy_volume(hdr,root_dir,format,metab_map,vSHIFT,outname)
% Output filename
%-------------------------------------------------------------------
fname = getfilelocation(hdr{1}, root_dir,'S',format);

% guess private field to use
if isfield(hdr{1}, 'Private_0029_1210')
    privdat = hdr{1}.Private_0029_1210;
elseif isfield(hdr{1}, 'Private_0029_1110')
    privdat = hdr{1}.Private_0029_1110;
else
    disp('Don''t know how to handle these spectroscopy data');
    fname = '';
    return;
end

% Image dimensions
%-------------------------------------------------------------------
nc = get_numaris4_numval(privdat,'Columns');
nr = get_numaris4_numval(privdat,'Rows');
% Guess number of timepoints in file - don't know whether this should be
% 'DataPointRows'-by-'DataPointColumns' or 'SpectroscopyAcquisitionDataColumns'
ntp = get_numaris4_numval(privdat,'DataPointRows')*get_numaris4_numval(privdat,'DataPointColumns');

%dim    = [nc nr numel(hdr) 2 ntp];
dim=[48 48 1];
nc=48;nr=48;
dt     = spm_type('float32'); % Fixed datatype

% Orientation information
%-------------------------------------------------------------------
% Axial Analyze voxel co-ordinate system:
% x increases     right to left
% y increases posterior to anterior
% z increases  inferior to superior

% DICOM patient co-ordinate system:
% x increases     right to left
% y increases  anterior to posterior
% z increases  inferior to superior

% T&T co-ordinate system:
% x increases      left to right
% y increases posterior to anterior
% z increases  inferior to superior

analyze_to_dicom = [diag([1 -1 1]) [0 (dim(2)+1) 0]'; 0 0 0 1]; % Flip voxels in y
patient_to_tal   = diag([-1 -1 1 1]); % Flip mm coords in x and y directions
shift_vx         = [eye(4,3) [.5; .5; 0; 1]];

orient           = reshape(get_numaris4_numval(privdat,...
                                               'ImageOrientationPatient'),[3 2]);
ps               = get_numaris4_numval(privdat,'PixelSpacing');
if nc*nr == 1
    % Single Voxel Spectroscopy (based on the following information from SIEMENS)
    %---------------------------------------------------------------
    % NOTE: Internally the position vector of the CSI matrix shows to the outer border
    % of the first voxel. Therefore the position vector has to be corrected.
    % (Note: The convention of Siemens spectroscopy raw data is in contrast to the
    %  DICOM standard where the position vector points to the center of the first voxel.)
    %---------------------------------------------------------------
    % SIEMENS decides which definition to use based on the contents of the
    % 'PixelSpacing' internal header field. If it has non-zero values,
    % assume DICOM convention. If any value is zero, assume SIEMENS
    % internal convention for this direction.
    % Note that in SIEMENS code, there is a shift when PixelSpacing is
    % zero. Here, the shift seems to be necessary when PixelSpacing is
    % non-zero. This may indicate more fundamental problems with
    % orientation decoding.
    if ps(1) == 0 % row
        ps(1) = get_numaris4_numval(privdat,...
                                    'VoiPhaseFoV');
        shift_vx(1,4) = 0;
    end
    if ps(2) == 0 % col
        ps(2) = get_numaris4_numval(privdat,...
                                    'VoiReadoutFoV');
        shift_vx(2,4) = 0;
    end
end
ps=[5;5];
pos = get_numaris4_numval(privdat,'ImagePositionPatient');
% for some reason, pixel spacing needs to be swapped
R  = [orient*diag(ps([2 1])); 0 0];
x1 = [1;1;1;1];
y1 = [pos; 1];

if length(hdr)>1,
    error('spm_dicom_convert:spectroscopy',...
        'Don''t know how to handle multislice spectroscopy data.');
else
    orient(:,3)      = null(orient');
    if det(orient)<0, orient(:,3) = -orient(:,3); end;
    try
        z = get_numaris4_numval(privdat,...
            'VoiThickness');
    catch
        try
            z = get_numaris4_numval(privdat,...
                'SliceThickness');
        catch
            z = 1;
        end
    end;
    x2 = [0;0;1;0];
    y2 = [orient*[0;0;z];0];
end
dicom_to_patient = [y1 y2 R]/[x1 x2 eye(4,2)];
mat              = patient_to_tal*dicom_to_patient*shift_vx*analyze_to_dicom;

% Possibly useful information
%-------------------------------------------------------------------
if checkfields(hdr{1},'AcquisitionTime','MagneticFieldStrength','MRAcquisitionType',...
        'ScanningSequence','RepetitionTime','EchoTime','FlipAngle',...
        'AcquisitionDate'),
    tim = datevec(hdr{1}.AcquisitionTime/(24*60*60));
    descrip = sprintf('%gT %s %s TR=%gms/TE=%gms/FA=%gdeg %s %d:%d:%.5g',...
        hdr{1}.MagneticFieldStrength, hdr{1}.MRAcquisitionType,...
        deblank(hdr{1}.ScanningSequence),...
        hdr{1}.RepetitionTime,hdr{1}.EchoTime,hdr{1}.FlipAngle,...
        datestr(hdr{1}.AcquisitionDate),tim(4),tim(5),tim(6));
else
    descrip = hdr{1}.Modality;
end;

if ~true, % LEFT-HANDED STORAGE
    mat    = mat*[-1 0 0 (dim(1)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
end;

% Write the image volume
%-------------------------------------------------------------------
N      = nifti;
pinfo  = [1 0];
if isfield(hdr{1},'RescaleSlope'),      pinfo(1) = hdr{1}.RescaleSlope;     end;
if isfield(hdr{1},'RescaleIntercept'),  pinfo(2) = hdr{1}.RescaleIntercept; end;

dim = [size(metab_map,1) size(metab_map,2) 1 size(metab_map,3)];
N.dat  = file_array(outname,dim,dt,0,pinfo(1),pinfo(2));
N.mat  = mat;
N.mat0 = mat;
%%
%%
%%
%% Here added
%%
%%
%%
cc = metab_map;
cc3 = imrotate(cc,-90);
%cc3 = fliplr(cc3);
 cc3 = circshift(cc3,vSHIFT(1),1);
 cc3 = circshift(cc3,vSHIFT(2),1);
% cc3 = circshift(cc3,-2,2);
%  cc3 = circshift(cc3,9,1); % sets ydir pixel shift.  Need to play to find best adjustment
temp1 = zeros(dim);
for i = 1:size(temp1,4)
    temp1(:,:,:,i,:) = cc3(:,:,i,:);
end
N.dat(:,:,:,:,:)=temp1;
%%
N.mat_intent  = 'Scanner';
N.mat0_intent = 'Scanner';
N.descrip     = descrip;
N.extras      = struct('MagneticFieldStrength',...
                       get_numaris4_numval(privdat,'MagneticFieldStrength'),...
                       'TransmitterReferenceAmplitude',...
                       get_numaris4_numval(privdat,'TransmitterReferenceAmplitude'));
create(N);

% Read data, swap dimensions
%data = permute(reshape(read_spect_data(hdr{1},privdat),dim([4 5 1 2 3])), ...
 %               [3 4 5 1 2]);
% plane = fliplr(plane);

%N.dat(:,:,:,:,:) = data;%uzay
return;
%_______________________________________________________________________

function [images,guff] = select_tomographic_images(hdr)
images = {};
guff   = {};
for i=1:length(hdr),
    if ~checkfields(hdr{i},'Modality') || ~(strcmp(hdr{i}.Modality,'MR') ||...
            strcmp(hdr{i}.Modality,'PT') || strcmp(hdr{i}.Modality,'CT'))
        if checkfields(hdr{i},'Modality'),
            fprintf('File "%s" can not be converted because it is of type "%s", which is not MRI, CT or PET.\n', hdr{i}.Filename, hdr{i}.Modality);
        else
            fprintf('File "%s" can not be converted because it does not encode an image.\n', hdr{i}.Filename);
        end
        guff = [guff(:)',hdr(i)];
    elseif ~checkfields(hdr{i},'StartOfPixelData','SamplesperPixel',...
            'Rows','Columns','BitsAllocated','BitsStored','HighBit','PixelRepresentation'),
        disp(['Cant find "Image Pixel" information for "' hdr{i}.Filename '".']);
        guff = [guff(:)',hdr(i)];
   %elseif isfield(hdr{i},'Private_2001_105f'),
   %    % This field corresponds to: > Stack Sequence 2001,105F SQ VNAP, COPY
   %    % http://www.medical.philips.com/main/company/connectivity/mri/index.html
   %    % No documentation about this private field is yet available.
   %    disp('Cant yet convert Phillips Intera DICOM.');
   %    guff = {guff{:},hdr{i}};
    elseif ~(checkfields(hdr{i},'PixelSpacing','ImagePositionPatient','ImageOrientationPatient')||isfield(hdr{i},'Private_0029_1110')||isfield(hdr{i},'Private_0029_1210')),
        disp(['Cant find "Image Plane" information for "' hdr{i}.Filename '".']);
        guff = [guff(:)',hdr(i)];
    elseif ~checkfields(hdr{i},'PatientID','SeriesNumber','AcquisitionNumber','InstanceNumber'),
       %disp(['Cant find suitable filename info for "' hdr{i}.Filename '".']);
        if ~isfield(hdr{i},'SeriesNumber')
            disp('Setting SeriesNumber to 1');
            hdr{i}.SeriesNumber=1;
            images = [images(:)',hdr(i)];
        end;
        if ~isfield(hdr{i},'AcquisitionNumber')
            if isfield(hdr{i},'Manufacturer') && ~isempty(strfind(upper(hdr{1}.Manufacturer), 'PHILIPS'))
                % WHY DO PHILIPS DO THINGS LIKE THIS????
                if isfield(hdr{i},'InstanceNumber')
                     hdr{i}.AcquisitionNumber = hdr{i}.InstanceNumber;
                else
                     disp('Setting AcquisitionNumber to 1');
                     hdr{i}.AcquisitionNumber=1;
                end
             else
                disp('Setting AcquisitionNumber to 1');
                hdr{i}.AcquisitionNumber=1;
             end
            images = [images(:)',hdr(i)];
        end;
        if ~isfield(hdr{i},'InstanceNumber')
            disp('Setting InstanceNumber to 1');
            hdr{i}.InstanceNumber=1;
            images = [images(:)',hdr(i)];
        end;
    else
        images = [images(:)',hdr(i)];
    end;
end;
return;
%%

function fname = getfilelocation(hdr,root_dir,prefix,format)

if nargin < 3
    prefix = 'f';
end;

if strncmp(root_dir,'ice',3)
    root_dir = root_dir(4:end);
    imtype = textscan(hdr.ImageType,'%s','delimiter','\\');
    try
        imtype = imtype{1}{3};
    catch
        imtype = '';
    end;
    prefix = [prefix imtype get_numaris4_val(hdr.CSAImageHeaderInfo,'ICE_Dims')];
end;

if strcmp(root_dir, 'flat')
    % Standard SPM file conversion
    %-------------------------------------------------------------------
    if checkfields(hdr,'SeriesNumber','AcquisitionNumber')
        if checkfields(hdr,'EchoNumbers')
            fname = sprintf('%s%s-%.4d-%.5d-%.6d-%.2d.%s', prefix, strip_unwanted(hdr.PatientID),...
                hdr.SeriesNumber, hdr.AcquisitionNumber, hdr.InstanceNumber,...
                hdr.EchoNumbers, format);
        else
            fname = sprintf('%s%s-%.4d-%.5d-%.6d.%s', prefix, strip_unwanted(hdr.PatientID),...
                hdr.SeriesNumber, hdr.AcquisitionNumber, ...
                hdr.InstanceNumber, format);
        end;
    else
        fname = sprintf('%s%s-%.6d.%s',prefix, ...
            strip_unwanted(hdr.PatientID),hdr.InstanceNumber, format);
    end;

    fname = fullfile(pwd,fname);
    return;
end;

% more fancy stuff - sort images into subdirectories
if ~isfield(hdr,'ProtocolName')
    if isfield(hdr,'SequenceName')
        hdr.ProtocolName = hdr.SequenceName;
    else
        hdr.ProtocolName='unknown';
    end;
end;
if ~isfield(hdr,'SeriesDescription')
    hdr.SeriesDescription = 'unknown';
end;
if ~isfield(hdr,'EchoNumbers')
    hdr.EchoNumbers = 0;
end;

m = sprintf('%02d', floor(rem(hdr.StudyTime/60,60)));
h = sprintf('%02d', floor(hdr.StudyTime/3600));
studydate = sprintf('%s_%s-%s', datestr(hdr.StudyDate,'yyyy-mm-dd'), ...
    h,m);
switch root_dir
    case {'date_time','series'}
    id = studydate;
    case {'patid', 'patid_date', 'patname'},
    id = strip_unwanted(hdr.PatientID);
end;
serdes = strrep(strip_unwanted(hdr.SeriesDescription),...
    strip_unwanted(hdr.ProtocolName),'');
protname = sprintf('%s%s_%.4d',strip_unwanted(hdr.ProtocolName), ...
    serdes, hdr.SeriesNumber);
switch root_dir
    case 'date_time',
        dname = fullfile(pwd, id, protname);
    case 'patid',
        dname = fullfile(pwd, id, protname);
    case 'patid_date',
        dname = fullfile(pwd, id, studydate, protname);
    case 'patname',
        dname = fullfile(pwd, strip_unwanted(hdr.PatientsName), ...
            id, protname);
    case 'series',
        dname = fullfile(pwd, protname);
    otherwise
        error('unknown file root specification');
end;
if ~exist(dname,'dir'),
    mkdir_rec(dname);
end;

% some non-product sequences on SIEMENS scanners seem to have problems
% with image numbering in MOSAICs - doublettes, unreliable ordering
% etc. To distinguish, always include Acquisition time in image name
sa = sprintf('%02d', floor(rem(hdr.AcquisitionTime,60)));
ma = sprintf('%02d', floor(rem(hdr.AcquisitionTime/60,60)));
ha = sprintf('%02d', floor(hdr.AcquisitionTime/3600));
fname = sprintf('%s%s-%s%s%s-%.5d-%.5d-%d.%s', prefix, id, ha, ma, sa, ...
        hdr.AcquisitionNumber,hdr.InstanceNumber, ...
        hdr.EchoNumbers,format);

fname = fullfile(dname, fname);

%_______________________________________________________________________
function ok = checkfields(hdr,varargin)
ok = 1;
for i=1:(nargin-1),
    if ~isfield(hdr,varargin{i}),
        ok = 0;
        break;
    end;
end;
return;

function [spect,images] = select_spectroscopy_images(hdr)
spectsel = zeros(1,numel(hdr));
for i=1:length(hdr),
    if isfield(hdr{i},'SOPClassUID')
        spectsel(i) = strcmp(hdr{i}.SOPClassUID,'1.3.12.2.1107.5.9.1');
    end;
end;
spect  = hdr(logical(spectsel));
images = hdr(~logical(spectsel));
return;
%___________

function clean = strip_unwanted(dirty)
msk = (dirty>='a'&dirty<='z') | (dirty>='A'&dirty<='Z') |...
      (dirty>='0'&dirty<='9') | dirty=='_';
clean = dirty(msk);
return;

function val = get_numaris4_numval(str,name)
val1 = get_numaris4_val(str,name);
val  = zeros(size(val1,1),1);
for k = 1:size(val1,1)
    val(k)=str2num(val1(k,:));
end;
return;

function val = get_numaris4_val(str,name)
name = deblank(name);
val  = {};
for i=1:length(str),
    if strcmp(deblank(str(i).name),name),
        for j=1:str(i).nitems,
            if  str(i).item(j).xx(1),
                val = {val{:} str(i).item(j).val};
            end;
        end;
        break;
    end;
end;
val = strvcat(val{:});
return;
%%