function [Ps,Pd,Pv] = FreemanDurden(C3)

arguments
    C3 (:,:,3,3)
end

fv = squeeze(4 * real(C3(:,:,2,2)));
A = squeeze(C3(:,:,1,3) - fv / 8);
B = squeeze(real(C3(:,:,3,3)) - 3 * fv / 8);
C = squeeze(real(C3(:,:,1,1)) - 3 * fv / 8);
D = conj(A);

% real(ShhSvv*) > 0
mask1 = real(squeeze(C3(:,:,1,3)) - fv / 8) > 0;
alpha1 = -1;
fd1 = real((B .* C - A .* D) ./ (A + B + C + D));
fs1 = B - fd1;
beta1 = (A + fd1) ./ fs1;

% real(ShhSvv*) <= 0
mask2 = real(squeeze(C3(:,:,1,3)) - fv / 8) <= 0;
beta2 = 1;
fs2 = real((A .* D - B .* C) ./ (A + D - B - C));
fd2 = B - fs2;
alpha2 = (A - fs2) ./ fd2;

fs = fs1 .* mask1 + fs2 .* mask2;
fd = fd1 .* mask1 + fd2 .* mask2;
alpha = alpha1 .* mask1 + alpha2 .* mask2;
beta = beta1 .* mask1 + beta2 .* mask2;

Ps = abs(fs .* (1 + abs(beta).^2));
Pd = abs(fd .* (1 + abs(alpha).^2));
Pv = abs(fv);

end