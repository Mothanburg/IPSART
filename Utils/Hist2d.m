% Hist2d - 2D 热力图
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function Hist2d(x, y, options)

arguments
    x (1,:)
    y (1,:)
    options.xLim (1,2) = [min(x) max(x)]
    options.yLim (1,2) = [min(y) max(y)]
    options.xNumBins = 200
    options.yNumBins = 200
    options.cmap string = "jet"
    options.background (1,3) = [1 1 1]
end

xbins = linspace(options.xLim(1), options.xLim(2), options.xNumBins);
ybins = linspace(options.yLim(1), options.yLim(2), options.yNumBins);

xi = round(interp1(xbins, 1:options.xNumBins, x, 'linear', 'extrap'));
yi = round(interp1(ybins, 1:options.yNumBins, y, 'linear', 'extrap'));
h = accumarray([yi(:),xi(:)], 1, [options.yNumBins options.xNumBins]);

imagesc(xbins, ybins, h); axis xy
cm = colormap(options.cmap);
cm(1,:) = options.background;
colormap(cm);

end