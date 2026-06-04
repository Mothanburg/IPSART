% CoherenceComputing - 干涉图相干性计算
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function coh = CoherenceComputing(master, slave, windowSize)

arguments
    master (:,:)
    slave (:,:)
    windowSize (1,1) {mustBeInteger} 
end

window = ones(windowSize) / windowSize^2;

nom = imfilter(master .* conj(slave), window);

filted = imfilter( ...
    abs(master).^2 + 1i * abs(slave).^2, ...
    window ...
    );
den = sqrt(real(filted) .* imag(filted));

coh = abs(nom) ./ den;

end