function result = RefinedLeeFilter(C3, nLooks)

arguments
    C3 PolC3
    nLooks {mustBeGreaterThanOrEqual(nLooks, 1)}
end

global GARS_CONFIG

if GARS_CONFIG.CAPABILITY >= 1
    try
        [errno, c11, c22, c33, c12r, c13r, c23r, c12i, c13i, c23i] = ...
            clib.gars.CloudePottier(C3.m11, C3.m22, C3.m33, C3.m12_r, C3.m13_r, ...
            C3.m23_r, C3.m12_i, C3.m13_i, C3.m23_i);
        if errno ~= 0
            error("error number %d returned.", errno);
        end
        result = PolC3(c11, c22, c33, c12r, c13r, c23r, c12i, c13i, c23i);
        return
    catch
        warning("An error occurred when calling library, fallback to matlab.\n" + ...
            "    Error message: '%s'", e.message());
    end
end

result = RLF_matlab(C3, nLooks);

end

function PolC3_out = RLF_matlab(PolC3_in, Looks)

dtype = PolC3_in.Dtype;
C3 = zeros(PolC3_in.Height, PolC3_in.Width, 3, 3, dtype);
for i = 1:3
    for j = 1:3
        C3(:,:,i, j) = PolC3_in.getPageAt(i, j);
    end
end

[height,width,dim,~] = size(C3);
span = squeeze(C3(:,:,1,1) + C3(:,:,2,2) + C3(:,:,3,3));
span_ex = parallel.pool.Constant([zeros(3, width+3*2, dtype);   % 在SPAN周围补一圈0
    zeros(height, 3, dtype), span, zeros(height, 3, dtype);
    zeros(3, width+3*2, dtype)]);

w1 = [-1,0,1;-1,0,1;-1,0,1];
w2 = [0,1,1;-1,0,1;-1,-1,0];
w3 = [1,1,1;0,0,0;-1,-1,-1];
w4 = [1,1,0;1,0,-1;0,-1,-1];
W = cat(3, w1, w2, w3, w4, fliplr(w1), w2', flipud(w3), flipud(w2));       % 边缘检测算子

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

Prewitt = parallel.pool.Constant({pw1, pw2, pw3, pw4, pw5, pw6, pw7, pw8});  % prewitt算子

C3 = parallel.pool.Constant(shiftdim(C3, 2));
result_shift = zeros(dim, dim, height, width, dtype);
parfor j = (1:width) + 3
    span_slice = span_ex.Value(:,(j - 3):(j + 3));
    result_col = zeros(dim, dim, height, dtype);
    for i = (1:height) + 3
        window = span_slice((i - 3):(i + 3),:);
        % 计算3x3的平均矩阵
        mx = zeros(3, 3, dtype);
        for p = 1:3
            for q = 1:3
                mx(p,q) = mean(window((1:3) + 2 * (p - 1),(1:3) + 2 * (q - 1)), "all");
            end
        end

        % 计算匹配的模板
        mx = mx .* W;
        [~,idx] = max(sum(mx, [1 2]));
        pw = Prewitt.Value{idx};

        % 计算b参数
        z = window .* pw;
        z_mean = sum(z, "all") / 28;
        var_v = 1 / Looks;
        var_z = sum((z - z_mean).^2, "all") / 28;
        var_x = (var_z - (z_mean^2) * var_v) / (1 + var_v);
        b = (var_x + 1e-30) / (var_z + 1e-30);

        % 计算协方差矩阵的统计平均值
        m = i - 3;
        n = j - 3;
        row = max([m-3, 1]):min([m+3, height]);
        if m < 4
            row_pw = 5-m:7;
        elseif m > height - 4
            row_pw = 1:4+(height-m);
        else
            row_pw = 1:7;
        end
        col = max([n-3, 1]):min([n+3, width]);
        if n < 4
            col_pw = 5-n:7;
        elseif n > width - 4
            col_pw = 1:4+(width-n);
        else
            col_pw = 1:7;
        end
        pw_masked = pw(row_pw,col_pw);
        c_mean = sum( ...
            C3.Value(:,:,row,col) .* reshape(pw_masked, [1 1 size(pw_masked)]), ...
            [3 4]) / sum(pw_masked, "all");

        c_caret = squeeze(c_mean + b * (C3.Value(:,:,i-3,j-3) - c_mean));
        result_col(:,:,i-3) = c_caret;
    end
    result_shift(:,:,:,j-3) = result_col;
end
result = shiftdim(result_shift, 2);

PolC3_out = PolC3(squeeze(result(:,:,1,1)), squeeze(result(:,:,2,2)), ...
    squeeze(result(:,:,3,3)), squeeze(real(result(:,:,1,2))), ...
    squeeze(real(result(:,:,1,3))), squeeze(real(result(:,:,2,3))), ...
    squeeze(imag(result(:,:,1,2))), squeeze(imag(result(:,:,1,3))), ...
    squeeze(imag(result(:,:,2,3))));

end
