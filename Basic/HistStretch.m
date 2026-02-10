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
    case "none"
        result = bands;
    case "linear"
        a = min(bands, [], [1 2]);
        b = max(bands, [], [1 2]);
        result = (bands - a) ./ (b - a);
    case "linear percent"
        narginchk(4, 4);
        a = prctile(bands, varargin{1}, [1 2]);
        b = prctile(bands, varargin{2}, [1 2]);
        result = (bands - a) ./ (b - a);
        result(result < 0) = 0;
        result(result > 1) = 1;
    case "optimized linear"
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
    case "log"
        result = log10(bands + 1e-45);
        a = min(result, [], [1 2]);
        b = max(result, [], [1 2]);
        result = (result - a) ./ (b - a);
    otherwise
        error("Unknown stretch method: %s", method);
end

result = squeeze(result);

end