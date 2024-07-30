function result = RefinedLeeFilter(M, nLook)

arguments
    M PolMat
    nLook {mustBeGreaterThanOrEqual(nLook, 1)}
end

global GARS_CONFIG

if GARS_CONFIG.CAPABILITY >= 1
    try
        if isa(M, "PolM3")
            [errno, m11, m22, m33, m12r, m13r, m23r, m12i, m13i, m23i] = ...
                clib.gars.RefinedLeeFilter3x3(nLook, M.m11, M.m22, M.m33, M.m12_r, M.m13_r, ...
                M.m23_r, M.m12_i, M.m13_i, M.m23_i);
            if errno ~= 0
                error("error number %d is returned.", errno);
            end

            if isa(M, "PolC3")
                result = PolC3(m11, m22, m33, m12r, m13r, m23r, m12i, m13i, m23i);
            elseif isa(M, "PolT3")
                result = PolT3(m11, m22, m33, m12r, m13r, m23r, m12i, m13i, m23i);
            end

            return;
        elseif isa(M, "PolM2")
            [errno, m11, m22, m12r, m12i] = ...
                clib.gars.RefinedLeeFilter2x2(nLook, M.m11, M.m22, M.m12_r, M.m12_i);
            if errno ~= 0
                error("error number %d is returned.", errno);
            end

            if isa(M, "PolC2")
                result = PolC2(m11, m22, m12r, m12i, M.PolType);
            elseif isa(M, "PolT2")
                result = PolT2(m11, m22, m12r, m12i);
            end

            return;
        end
    catch e
        warning(e.identifier, "An error occurred when calling library, fallback to matlab.\n" + ...
            "        Error message: %s", e.message);
    end
end

result = RLF_matlab(M, nLook);

end

function result = RLF_matlab(M, nLook)

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

height = M.Height;
width = M.Width;

pw_id = zeros(height, width, "int32");
b = zeros(height, width, M.Dtype);

span = parallel.pool.Constant(M.SPAN);
M = parallel.pool.Constant(M);
parfor col = 1:width
    for row = 1:height
        window = zeros(7, 7, M.Value.Dtype);
        for i = 1:7
            for j = 1:7
                r = min(max(row - 3 + i, 1), height);
                c = min(max(col - 3 + j, 1), width);
                window(i,j) = span.Value(r, c);
            end
        end

        meanmat = zeros(3, 3, M.Value.Dtype);
        for i = 1:3
            for j = 1:3
                meanmat(i,j) = mean(window((1:3) + 2 * (i - 1),(1:3) + 2 * (j - 1)), "all");
            end
        end
        [~,wid] = max(sum(meanmat .* W.Value, [1 2]));
        pw_id(row,col) = wid;

        pw = PW.Value{wid};
        z = window .* pw;
        z_mean = sum(z, "all") / 28;
        var_v = 1 / nLook;
        var_z = sum((z - z_mean).^2, "all") / 28;
        var_x = (var_z - (z_mean^2) * var_v) / (1 + var_v);
        b(row,col) = (var_x + 1e-30) / (var_z + 1e-30);
    end
end

result = M.Value.fmapPage(UnaryOp(@rlf_filt, 1, pw_id, b));

end

function out = rlf_filt(cij, pwId, b)

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

[height,width] = size(cij);
out = zeros(height, width, class(cij));

cij = parallel.pool.Constant(cij);
b = parallel.pool.Constant(b);
parfor col = 1:width
    for row = 1:height
        pw = PW.Value{pwId};

        m_mean = cast(0, "like", cij.Value);
        for i = 1:7
            for j = 1:7
                r = min(max(row - 3 + i, 1), height);
                c = min(max(col - 3 + j, 1), width);
                m_mean = m_mean + cij.Value(r, c) * pw(i, j);
            end
        end
        m_mean = m_mean / 28;
        m_out = m_mean + b.Value(row, col) * (cij.Value(row, col) - m_mean);

        out(row, col) = m_out;
    end
end


end

