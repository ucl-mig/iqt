function dt= BigbassReadDT_Volume(nameroot)

for i=1:8
    dt(:,:,:,i) = read_std_nii([ nameroot num2str(i) '.nii' ]);
end
