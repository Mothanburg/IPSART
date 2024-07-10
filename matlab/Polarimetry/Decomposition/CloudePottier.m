% Cloude-Pottier H/a/A decomposition
function [H,alpha,A] = CloudePottier(M3, options)

arguments
    M3 PolM3
    options.quiet logical = false
end

if ~options.quiet && ~isa(M3, "PolT3")
    warning("The input data is not the coherency matrix T, the returned alpha may not make sense");
end

global GARS_CONFIG

if GARS_CONFIG.CAPABILITY > 0
    try
        [H,alpha,A] = clib.ga9rs.CloudePottier(M3.m11, M3.m22, M3.m33, M3.m12_r, ...
            M3.m13_r, M3.m23_r, M3.m12_i, M3.m13_i, M3.m23_i);
        return
    catch e
        warning(e.identifier, "An error occurred when calling library, fallback to matlab.\n" + ...
            "        Error message: %s", e.message);
    end
end

[H,alpha,A] = CP_matlab(M3);

end


function [H,alpha,A] = CP_matlab(M3)

height = M3.Height;
width = M3.Width;

H = zeros(height, width, M3.Dtype);
alpha = zeros(height, width, M3.Dtype);
A = zeros(height, width, M3.Dtype);

M3 = parallel.pool.Constant(M3);
parfor j = 1:width
    for i = 1:height
        t = M3.Value.getMatAt(i, j);
        [v,d] = eig(t);
        d = diag(abs(d))';
        p = d / sum(d);
        H(i,j) = calc_entropy(p);
        alpha(i,j) = sum(p .* acosd(abs(v(1,:))));
        l3 = min(d);
        l2 = sum(d) - max(d) - l3;
        A(i,j) = (l2 - l3) / (l2 + l3);
    end
end

end

function h = calc_entropy(x)

h = 0;
len = numel(x);
for i = 1:len
    if x ~= 0
        h = h - x(i) * log(x(i)) / log(len);
    end
end

end