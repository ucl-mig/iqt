function [dt, hdr] = ReadDT_Volume(nameroot)
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%
[~, hdr] = read_std_nii([ nameroot '1.nii' ]);
for i=1:8
    dt(:,:,:,i) = read_std_nii([ nameroot num2str(i) '.nii' ]);
end
