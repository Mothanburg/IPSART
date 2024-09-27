function [H,alpha,A] = CloudePottierDP(M2)

arguments
    M2 PolM2
end

global gars_cp2_cpu_enable
if isempty(gars_cp2_cpu_enable)
    gars_cp2_cpu_enable = true;
end

if gars_cp2_cpu_enable
    try
        [H,alpha,A] = CP2_native(M2);
        return;
    catch e
        warning(e.identifier, ...
            "An error occurred when calling library, fallback to matlab\n" + ...
            "        Error message: %s", e.message);
        gars_cp2_cpu_enable = false;
    end
end

[H,alpha,A] = CP2_matlab(M2);

end


%---------- native function caller ----------%

function [H,alpha,A] = CP2_native(M2)

if M2.Dtype == "double"
    [H,alpha,A] = clib.gars.CloudePottier2d(M2.m11, M2.m22, M2.m12_r, M2.m12_i);
else
    [H,alpha,A] = clib.gars.CloudePottier2f(M2.m11, M2.m22, M2.m12_r, M2.m12_i);
end

end


%---------- MATLAB function caller ----------%

function [H,alpha,A] = CP2_matlab(M2)

height = M2.Height;
width = M2.Width;

H = zeros(height, width, M2.Dtype);
alpha = zeros(height, width, M2.Dtype);
A = zeros(height, width, M2.Dtype);

M2 = parallel.pool.Constant(M2);
parfor j = 1:width
    for i = 1:height
        t = M2.Value.MatAt(i, j);
        [v,d] = eig(t);
        d = diag(abs(d))';
        p = d / sum(d);
        H(i,j) = calc_entropy(p);
        alpha(i,j) = sum(p .* acosd(abs(v(1,:))));
        A(i,j) = abs(d(1) - d(2)) / sum(d);
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