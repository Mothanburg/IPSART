function C3 = fn_gen_C3(HH, HV, VH, VV)

X = (HV + VH) / sqrt(2);
C3(:,:,1,1) = HH .* conj(HH);
C3(:,:,1,2) = HH .* conj(X);
C3(:,:,1,3) = HH .* conj(VV);
C3(:,:,2,1) = conj(C3(:,:,1,2));
C3(:,:,2,2) = X .* conj(X);
C3(:,:,2,3) = X .* conj(VV);
C3(:,:,3,1) = conj(C3(:,:,1,3));
C3(:,:,3,2) = conj(C3(:,:,2,3));
C3(:,:,3,3) = VV .* conj(VV);

end