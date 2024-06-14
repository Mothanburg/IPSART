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
if length(sz) < 3
    image = reshape(image, [sz, 1]);
    sz = [sz 1];
end

in_rows = sz(1);
in_cols = sz(2);

row_strides = repmat(rowLook, 1, floor(in_rows / rowLook));
row_left = rem(in_rows, rowLook);
if row_left == 0
    row_left = [];
end

col_strides = repmat(colLook, 1, floor(in_cols / colLook));
col_left = rem(in_cols, colLook);
if col_left == 0
    col_left = [];
end

image = parallel.pool.Constant(image);
dim_cnt = prod(sz(3:end));
parfor idx = 1:dim_cnt
    patched = mat2cell(...
        squeeze(image.Value(:,:,idx)), ...
        [row_strides row_left], ...
        [col_strides col_left] ...
        );
    result(:,:,idx) = cellfun(@(x) mean(x(:)), patched);
end

out_sz = size(result);
result = squeeze(reshape(result, [out_sz(1:2) sz(3:end)]));

end
