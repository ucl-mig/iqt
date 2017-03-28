function [I,elsp,Meta] = read_std_any(filename)

convert_image(filename,'tmp.nii');

%I = load_nii(filename);
I = load_untouch_nii('tmp.nii');

delete('tmp.nii');

Meta = I.hdr;
I = I.img;

elsp = Meta.dime.pixdim(2:4);

%I = flipdim(flipdim(I.img,1),2);
