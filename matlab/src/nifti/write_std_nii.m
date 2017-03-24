function write_std_nii(D,filename,spacing)

if (nargin<3)
    spacing = [1 1 1];
end

nii = make_nii(D);
nii.hdr.hist.qform_code = 0;
nii.hdr.hist.sform_code = 0;
nii.hdr.hist.magic = 'n+1';
nii.untouch = 1;
nii.hdr.dime.pixdim(2:4) = spacing;
save_untouch_nii(nii,filename);
