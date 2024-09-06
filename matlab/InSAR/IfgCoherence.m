function coh = IfgCoherence(master, slave, windowSize)

window = ones(windowSize);

nom = imfilter(master .* conj(slave), window);

filted = imfilter( ...
    master .* conj(master) + 1i * slave .* conj(slave), ...
    window ...
    );
den = sqrt(real(filted) .* imag(filted));

coh = abs(nom ./ den);

end