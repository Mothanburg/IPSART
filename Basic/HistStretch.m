% HistStretch - SAR图像亮度量化拉伸 
% 量化缩放SAR图像的像素亮度，使其便于显示
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function result = HistStretch(bands, method, varargin)

arguments
    bands (:,:,:) {mustBeReal}
    method string
end

arguments (Repeating)
    varargin {mustBeNumeric}
end

if length(size(bands)) < 3
    bands = reshape(bands, [size(bands) 1]);
end

switch lower(method)
    case "none" % 不做任何处理
        result = bands;
    case "linear" % 仅将像素值缩放到[0,1]之间
        a = min(bands, [], [1 2]);
        b = max(bands, [], [1 2]);
        result = (bands - a) ./ (b - a);
    case "linear percent" % 指定纯黑和纯白的像素亮度百分比阈值，并将像素值缩放到[0,1]之间
        narginchk(4, 4);
        a = prctile(bands, varargin{1}, [1 2]);
        b = prctile(bands, varargin{2}, [1 2]);
        result = (bands - a) ./ (b - a);
        result(result < 0) = 0;
        result(result > 1) = 1;
    case "optimized linear" % 优化的线性亮度缩放，并将像素值缩放到[0,1]之间
        a = prctile(bands, 2.5, [1 2]);
        b = prctile(bands, 99, [1 2]);
        c = a - 0.1 * (b - a);
        d = b + 0.5 * (b - a);
        [~,~,n_bands] = size(bands);
        for idx = 1:n_bands
            tmp = squeeze(bands(:,:,idx));
            tmp(tmp < c(idx)) = c(idx);
            tmp(tmp > d(idx)) = d(idx);
            result(:,:,idx) = tmp;
        end
        result = (result - c) ./ (d - c);
    case "log" % 对像素亮度取对数，并将像素值缩放到[0,1]之间
        result = log10(bands + 1e-45);
        a = min(result, [], [1 2]);
        b = max(result, [], [1 2]);
        result = (result - a) ./ (b - a);
    otherwise
        error("Unknown stretch method: %s", method);
end

result = squeeze(result);

end