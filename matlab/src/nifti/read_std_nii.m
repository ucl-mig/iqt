function [I,Meta] = read_std_nii(filename)

%I = load_nii(filename);
I = load_untouch_nii(filename);
    
Meta = I.hdr;
I = I.img;

%I = flipdim(flipdim(I.img,1),2);
