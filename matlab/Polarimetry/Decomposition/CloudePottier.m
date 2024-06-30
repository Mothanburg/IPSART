function [H,alpha,A] = CloudePottier(T3)

arguments
    T3 PolT3
end

global GARS_CONFIG

if GARS_CONFIG.CAPABILITY > 0
    try
        [H,alpha,A] = clib.gars.CloudePottier(T3.m11, T3.m22, T3.m33, T3.m12_r, ...
            T3.m13_r, T3.m23_r, T3.m12_i, T3.m13_i, T3.m23_i);
        return
    catch
        warning("An error occurred when calling library, fallback to matlab.\n" + ...
            "    Error message: '%s'", e.message());
    end
end

[H,alpha,A] = CP_matlab(T3);

end


function [H,alpha,A] = CP_matlab(T3)

height = T3.Height;
width = T3.Width;

H = zeros(height, width, T3.Dtype);
alpha = zeros(height, width, T3.Dtype);
A = zeros(height, width, T3.Dtype);

T3 = parallel.pool.Constant(T3);
parfor j = 1:width
    for i = 1:height
        t = T3.Value.getMatAt(i, j);
        [v,d] = eig(t);
        d = diag(abs(d))'
        p = d / sum(d);
        H(i,j) = -sum(p .* log(p) / log(3));
        alpha(i,j) = sum(p .* acosd(abs(v(1,:))));
        l3 = min(d);
        l2 = sum(d) - max(d) - l3;
        A(i,j) = (l2 - l3) / (l2 + l3);
    end
end

end