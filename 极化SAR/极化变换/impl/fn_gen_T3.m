function T3 = fn_gen_T3(HH, HV, VH, VV)

X = (HV + VH) / 2;
T3(:,:,1,1) = (HH + VV) .* conj(HH + VV);
T3(:,:,1,2) = (HH + VV) .* conj(HH - VV);
T3(:,:,1,3) = 2 * (HH + VV) .* conj(X);
T3(:,:,2,1) = conj(T3(:,:,1,2));
T3(:,:,2,2) = (HH - VV) .* conj(HH - VV);
T3(:,:,2,3) = 2 * (HH - VV) .* conj(X);
T3(:,:,3,1) = conj(T3(:,:,1,3));
T3(:,:,3,2) = conj(T3(:,:,2,3));
T3(:,:,3,3) = 4 * X .* conj(X);
T3 = T3 / 2;

end