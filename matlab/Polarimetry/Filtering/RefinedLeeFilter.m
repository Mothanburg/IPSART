function result = RefinedLeeFilter(M3, nLooks)

arguments
    M3 PolM3
    nLooks {mustBeGreaterThanOrEqual(nLooks, 1)}
end

global GARS_CONFIG

fallback_flag = false;
if GARS_CONFIG.CAPABILITY >= 1
    try
        [errno, m11, m22, m33, m12r, m13r, m23r, m12i, m13i, m23i] = ...
            clib.gars.RefinedLeeFilter(nLooks, M3.m11, M3.m22, M3.m33, M3.m12_r, M3.m13_r, ...
            M3.m23_r, M3.m12_i, M3.m13_i, M3.m23_i);
        if errno ~= 0
            error("error number %d is returned.", errno);
        end
    catch e
        warning(e.identifier, "An error occurred when calling library, fallback to matlab.\n" + ...
            "    Error message: %s", e.message);
        fallback_flag = true;
    end
end

if fallback_flag
    [m11, m22, m33, m12r, m13r, m23r, m12i, m13i, m23i] = RLF_matlab(M3, nLooks);
end

if isa(M3, "PolC3")
    result = PolC3(m11, m22, m33, m12r, m13r, m23r, m12i, m13i, m23i);
else
    result = PolT3(m11, m22, m33, m12r, m13r, m23r, m12i, m13i, m23i);
end

end

function [m11, m22, m33, m12r, m13r, m23r, m12i, m13i, m23i] = RLF_matlab(M3, nLook)

w1 = [-1,0,1;-1,0,1;-1,0,1];
w2 = [0,1,1;-1,0,1;-1,-1,0];
w3 = [1,1,1;0,0,0;-1,-1,-1];
w4 = [1,1,0;1,0,-1;0,-1,-1];
W = parallel.pool.Constant(cat(3, w1, w2, w3, w4, fliplr(w1), w2', flipud(w3), flipud(w2)));

pw1 = repmat([0,0,0,1,1,1,1], [7,1]);
pw2 = [1,1,1,1,1,1,1;
    0,1,1,1,1,1,1;
    0,0,1,1,1,1,1;
    0,0,0,1,1,1,1;
    0,0,0,0,1,1,1;
    0,0,0,0,0,1,1;
    0,0,0,0,0,0,1];
pw3 = repmat([1,1,1,1,0,0,0]', [1,7]);
pw4 = flip(pw2, 2);
pw5 = flip(pw1, 2);
pw6 = flip(pw4, 1);
pw7 = flip(pw3, 1);
pw8 = flip(pw2, 1);
PW = parallel.pool.Constant({pw1, pw2, pw3, pw4, pw5, pw6, pw7, pw8});

height = M3.Height;
width = M3.Width;

m11 = zeros(height, width, M3.Dtype);
m22 = zeros(height, width, M3.Dtype);
m33 = zeros(height, width, M3.Dtype);
m12r = zeros(height, width, M3.Dtype);
m13r = zeros(height, width, M3.Dtype);
m23r = zeros(height, width, M3.Dtype);
m12i = zeros(height, width, M3.Dtype);
m13i = zeros(height, width, M3.Dtype);
m23i = zeros(height, width, M3.Dtype);

span = parallel.pool.Constant(M3.SPAN);
M3 = parallel.pool.Constant(M3);
parfor col = 1:width
    for row = 1:height
        window = zeros(7, 7, M3.Value.Dtype);
        for i = 1:7
            for j = 1:7
                r = min(max(row - 3 + i, 1), height);
                c = min(max(col - 3 + j, 1), width);
                window(i,j) = span.Value(r, c);
            end
        end

        meanmat = zeros(3, 3, M3.Value.Dtype);
        for i = 1:3
            for j = 1:3
                meanmat(i,j) = mean(window((1:3) + 2 * (i - 1),(1:3) + 2 * (j - 1)), "all");
            end
        end
        [~,wid] = max(sum(meanmat .* W.Value, [1 2]));
        pw = PW.Value{wid};

        % 计算b参数
        z = window .* pw;
        z_mean = sum(z, "all") / 28;
        var_v = 1 / nLook;
        var_z = sum((z - z_mean).^2, "all") / 28;
        var_x = (var_z - (z_mean^2) * var_v) / (1 + var_v);
        b = (var_x + 1e-30) / (var_z + 1e-30);

        % 滤波
        m_mean = zeros(3, 3, M3.Value.Dtype);
        for i = 1:7
            for j = 1:7
                r = min(max(row - 3 + i, 1), height);
                c = min(max(col - 3 + j, 1), width);
                m_mean = m_mean + M3.Value.getMatAt(r, c) * pw(i, j);
            end
        end
        m_mean = m_mean / 28;
        m_out = m_mean + b * (M3.Value.getMatAt(row, col) - m_mean);

        m11(row, col) = real(m_out(1, 1));
        m22(row, col) = real(m_out(2, 2));
        m33(row, col) = real(m_out(3, 3));
        m12r(row, col) = real(m_out(1, 2));
        m13r(row, col) = real(m_out(1, 3));
        m23r(row, col) = real(m_out(2, 3));
        m12i(row, col) = imag(m_out(1, 2));
        m13i(row, col) = imag(m_out(1, 3));
        m23i(row, col) = imag(m_out(2, 3));
    end
end

end

