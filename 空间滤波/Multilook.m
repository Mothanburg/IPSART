% 多视
function result = Multilook(image, rowLook, colLook)

arguments
    image double
    rowLook double
    colLook double
end

sz = size(image);
if length(sz) < 2
    error("输入值至少必须是二维图像");
end

in_rows = sz(1);
in_cols = sz(2);

row_strides = 1:rowLook:in_rows;
col_strides = 1:colLook:in_cols;

out_rows = length(row_strides);
out_cols = length(col_strides);

if length(sz) < 3
    image = reshape(image, [sz, 1]);
    result = zeros(out_rows, out_cols, 1);
else
    result = zeros([out_rows out_cols sz(3:end)]);
end

row_strides = parallel.pool.Constant(row_strides);
col_strides = parallel.pool.Constant(col_strides);
image = parallel.pool.Constant(image);
parfor j = 1:out_cols
    for i = 1:out_rows
        row1 = row_strides.Value(i);
        if i == out_rows
            row2 = in_rows;
        else
            row2 = row_strides.Value(i + 1) - 1;
        end
        col1 = col_strides.Value(j);
        if j == out_cols
            col2 = in_cols;
        else
            col2 = col_strides.Value(j + 1) - 1;
        end
        result(i,j,:) = mean(image.Value(row1:row2,col1:col2,:), [1 2]);
    end
end

result = squeeze(result);

end
