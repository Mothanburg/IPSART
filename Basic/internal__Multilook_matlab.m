function result = internal__Multilook_matlab( ...
    image, in_rows, in_cols, row_look, col_look)

row_strides = repmat(row_look, 1, floor(in_rows / row_look));
row_left = rem(in_rows, row_look);

col_strides = repmat(col_look, 1, floor(in_cols / col_look));
col_left = rem(in_cols, col_look);

patched = mat2cell(image(1:end-row_left,1:end-col_left), ...
    row_strides, col_strides);

result = cellfun(@(x) mean(x(:)), patched);

end
