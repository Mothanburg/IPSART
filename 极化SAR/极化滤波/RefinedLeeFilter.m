% 精致Lee滤波的极化降斑
function result = RefinedLeeFilter(C, Looks)

arguments
    C (:,:,:,:) double
    Looks double {mustBeGreaterThanOrEqual(Looks, 1)}
end

[height,width,dim,~] = size(C);
span = zeros(height, width);
for idx = 1:dim
    span = span + squeeze(abs(C(:,:,idx,idx)));
end
span_ex = parallel.pool.Constant([ ...
    span(1,1) span(1,:) span(end,end);
    span(:,1) span(:,:) span(:,end);
    span(end,1) span(end,:) span(end,end)...
    ]); % 在span周围镜像填充一圈

w1 = [-1,0,1;-1,0,1;-1,0,1];
w2 = [0,1,1;-1,0,1;-1,-1,0];
w3 = [1,1,1;0,0,0;-1,-1,-1];
w4 = [1,1,0;1,0,-1;0,-1,-1];
W = cat(3, w1, w2, w3, w4);       % 边缘检测算子

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

C = parallel.pool.Constant(shiftdim(C, 2));
result_shift = zeros(dim, dim, height, width);
parfor j = (1:width) + 3
    span_slice = span_ex.Value(:,(j - 3):(j + 3));
    result_col = zeros(dim, dim, height);
    for i = (1:height) + 3
        window = span_slice((i - 3):(i + 3),:);
        % 计算3x3的平均矩阵
        m = zeros(3, 3);
        for p = 1:3
            for q = 1:3
                m(p,q) = mean(window((1:3) + 2 * (i - 1),(1:3) + 2 * (j - 1)), "all");
            end
        end

        % 计算匹配的模板
        m = m .* W;
        [~,idx] = max(abs(sum(m, [1 2])));
        switch idx
            case 1
                delta1 = abs(m(2,1) - m(2,2));
                delta2 = abs(m(2,3) - m(2,2));
                if delta1 > delta2
                    idx = idx + 4;
                end
            case 2
                delta1 = abs(m(3,1) - m(2,2));
                delta2 = abs(m(1,3) - m(2,2));
                if delta1 > delta2
                    idx = idx + 4;
                end
            case 3
                delta1 = abs(m(3,2) - m(2,2));
                delta2 = abs(m(1,2) - m(2,2));
                if delta1 > delta2
                    idx = idx + 4;
                end
            case 4
                delta1 = abs(m(3,3) - m(2,2));
                delta2 = abs(m(1,1) - m(2,2));
                if delta1 > delta2
                    idx = idx + 4;
                end
        end
        pw = Prewitt.Value{idx};

        % 计算b参数
        z = window .* pw;
        z_mean = mean(z, "all");
        var_v = 1 / Looks;
        var_z = mean((z - z_mean).^2, "all");
        var_x = (var_z - (z_mean^2) * var_v) / (1 + var_v);
        b = var_x / var_z;
        if (isnan(b))
            b = 0;
        end

        % 计算协方差矩阵的统计平均值
        row = max([i - 6, 1]):min([i, height]);
        if i < 7
            row_pw = (2 - i):7;
        elseif i > height - 1
            row_pw = 1:(7 + height - i);
        else
            row_pw = 1:7;
        end
        col = max([j - 6, 1]):min([j, width]);
        if j < 7
            col_pw = (2 - j):7;
        elseif j > width - 1
            col_pw = 1:(7 + width - j);
        else
            col_pw = 1:7;
        end
        pw_masked = pw(row_pw,col_pw);
        c_mean = mean(C.Value(:,:,row,col) .* reshape(pw_masked, [1 1 size(pw_masked)]), [3 4]);

        c_caret = squeeze(c_mean + b * (C.Value(:,:,i-3,j-3) - c_mean));
        result_col(:,:,i-3) = c_caret;
    end
    result_shift(:,:,:,j-3) = result_col;
end
result = shiftdim(result_shift, 2);

end
