function result = RefinedLeeFilter(C3, nLooks)

arguments
    C3 PolC3
    nLooks {mustBeGreaterThanOrEqual(nLooks, 1)}
end

global GARS_CONFIG

if GARS_CONFIG.CAPABILITY >= 1
    try
        [errno, c11, c22, c33, c12r, c13r, c23r, c12i, c13i, c23i] = ...
            clib.gars.RefinedLeeFilter(nLooks, C3.m11, C3.m22, C3.m33, C3.m12_r, C3.m13_r, ...
            C3.m23_r, C3.m12_i, C3.m13_i, C3.m23_i);
        if errno ~= 0
            error("error number %d is returned.", errno);
        end
        result = PolC3(c11, c22, c33, c12r, c13r, c23r, c12i, c13i, c23i);
        return
    catch e
        warning(e.identifier, "An error occurred when calling library, fallback to matlab.\n" + ...
            "    Error message: %s", e.message);
    end
end

result = RLF_matlab(C3, nLooks);

end

function result = RLF_matlab(C3, nLook)

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

height = C3.Height;
width = C3.Width;

c11 = zeros(height, width, C3.Dtype);
c22 = zeros(height, width, C3.Dtype);
c33 = zeros(height, width, C3.Dtype);
c12r = zeros(height, width, C3.Dtype);
c13r = zeros(height, width, C3.Dtype);
c23r = zeros(height, width, C3.Dtype);
c12i = zeros(height, width, C3.Dtype);
c13i = zeros(height, width, C3.Dtype);
c23i = zeros(height, width, C3.Dtype);

span = parallel.pool.Constant(C3.SPAN);
C3 = parallel.pool.Constant(C3);
parfor col = 1:width
    for row = 1:height
        window = zeros(7, 7, C3.Value.Dtype);
        for i = 1:7
            for j = 1:7
                r = min(max(row - 3 + i, 1), height);
                c = min(max(col - 3 + j, 1), width);
                window(i,j) = span.Value(r, c);
            end
        end

        meanmat = zeros(3, 3, C3.Value.Dtype);
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
        c_mean = zeros(3, 3, C3.Value.Dtype);
        for i = 1:7
            for j = 1:7
                r = min(max(row - 3 + i, 1), height);
                c = min(max(col - 3 + j, 1), width);
                c_mean = c_mean + C3.Value.getMatAt(r, c) * pw(i, j);
            end
        end
        c_mean = c_mean / 28;
        c_out = c_mean + b * (C3.Value.getMatAt(row, col) - c_mean);

        c11(row, col) = real(c_out(1, 1));
        c22(row, col) = real(c_out(2, 2));
        c33(row, col) = real(c_out(3, 3));
        c12r(row, col) = real(c_out(1, 2));
        c13r(row, col) = real(c_out(1, 3));
        c23r(row, col) = real(c_out(2, 3));
        c12i(row, col) = imag(c_out(1, 2));
        c13i(row, col) = imag(c_out(1, 3));
        c23i(row, col) = imag(c_out(2, 3));
    end
end

result = PolC3(c11, c22, c33, c12r, c13r, c23r, c12i, c13i, c23i);

end
