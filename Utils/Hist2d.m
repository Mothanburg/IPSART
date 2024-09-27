function Hist2d(x, y, xlim, ylim, xNumBins, yNumBins, options)

arguments
    x (1,:)
    y (1,:)
    xlim (1,2) = [min(x) max(x)]
    ylim (1,2) = [min(y) max(y)]
    xNumBins = 200
    yNumBins = 200
    options.cmap string = "jet"
    options.background (1,3) = [1 1 1]
end

xbins = linspace(xlim(1), xlim(2), xNumBins);
ybins = linspace(ylim(1), ylim(2), yNumBins);

xi = round(interp1(xbins, 1:xNumBins, x, 'linear', 'extrap'));
yi = round(interp1(ybins, 1:yNumBins, y, 'linear', 'extrap'));
h = accumarray([yi(:),xi(:)], 1, [yNumBins xNumBins]);

imagesc(xbins, ybins, h); axis xy
cm = colormap(options.cmap);
cm(1,:) = options.background;
colormap(cm);

end