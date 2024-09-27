function result = HistStretch(bands, method, varargin)

arguments
    bands (:,:,:) {mustBeReal}
    method string
end

arguments(Repeating)
    varargin {mustBeNumeric}
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
        [~,~,n_bands] = size(bands);
        for idx = 1:n_bands
            tmp = squeeze(bands(:,:,idx));
            tmp(tmp < c(idx)) = c(idx);
            tmp(tmp > d(idx)) = d(idx);
            result(:,:,idx) = tmp;
        end
        result = (result - c) ./ (d - c);
    case "Log"
        a = min(bands, [], [1 2]);
        result = 10 * log10(bands - a + 1);
        result = result ./ max(result, [], [1 2]);
    otherwise
        error("Unknown stretch method: %s", method);
end

result = squeeze(result);

end