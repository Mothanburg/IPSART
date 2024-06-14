function result = HistStretch(bands, method, varargin)

arguments
    bands (:,:,:) double {mustBeReal}
    method string
end

arguments(Repeating)
    varargin double
end

if length(size(bands)) < 3
    bands = reshape(bands, [size(bands) 1]);
end

switch method
    case "None"
        result = bands;
    case "Linear"
        a = min(bands, [], [1 2]);
        b = max(bands, [], [1 2]);
        result = (bands - a) ./ (b - a);
    case "Linear Percent"
        narginchk(4, 4);
        a = prctile(bands, varargin{1}, [1 2]);
        b = prctile(bands, varargin{2}, [1 2]);
        result = (bands - a) ./ (b - a);
        result(result < 0) = 0;
        result(result > 1) = 1;
    case "Optimized Linear"
        a = prctile(bands, 2.5, [1 2]);
        b = prctile(bands, 99, [1 2]);
        c = a - 0.1 * (b - a);
        d = b + 0.5 * (b - a);
        result = bands;
        result(bands < c) = c;
        result(bands > d) = d;
        result = (result - c) ./ (d - c);
    case "Log"
        a = min(bands, [], [1 2]);
        result = log10(bands - a + 1);
        result = result ./ max(result, [], [1 2]);
    otherwise
        error("HistStretch:unknownStretchMethod", "未知的拉伸方法：%s", method);
end

result = squeeze(result);

end