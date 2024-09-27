function G = c2_2_sv(C2)

[height,width,~,~] = size(C2);
G = zeros(height, width, 4);
G(:,:,1) = C2(:,:,1,1)+C2(:,:,2,2);
G(:,:,2) = C2(:,:,1,1)-C2(:,:,2,2);
G(:,:,3) = 2*real(C2(:,:,1,2));
G(:,:,4) = -2*imag(C2(:,:,1,2));

end