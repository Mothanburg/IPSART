function coh = IfgCoherence(master, slave, windowSize)

window = ones(windowSize);

nom = imfilter(master .* conj(slave), window);

filted = imfilter( ...
    abs(master).^2 + 1i * abs(slave).^2, ...
    window ...
    );
den = sqrt(real(filted) .* imag(filted));

coh = abs(nom) ./ den;

end