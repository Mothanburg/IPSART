% 多视
function result = Multilook(image, rowLook, colLook)

arguments
    image (:,:)
    rowLook {mustBeInteger}
    colLook {mustBeInteger}
end

sz = size(image);
if length(sz) < 2
    error("输入值至少必须是二维图像");
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

patched = mat2cell(...
    image, ...
    [row_strides row_left], ...
    [col_strides col_left] ...
    );
result(:,:) = cellfun(@(x) mean(x(:)), patched);

end
