function write_hdr_nii(D,filename,hdr)

nii = make_nii(D);
nii.hdr.hist.qform_code = 0;
nii.hdr.hist.sform_code = 0;
nii.hdr.hist.magic = 'n+1';
nii.untouch = 1;

if (nargin>2)
    hdr.dime.dim = nii.hdr.dime.dim;
    nii.hdr = hdr;
end

save_untouch_nii(nii,filename);
