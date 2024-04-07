function S2 = fn_gen_S2(HH, HV, VH, VV)

S2(:,:,1,1) = HH;
S2(:,:,1,2) = HV;
S2(:,:,2,1) = VH;
S2(:,:,2,2) = VV;

end