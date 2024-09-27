function result = Multilook(image, rowLook, colLook)

arguments
    image (:,:)
    rowLook {mustBeInteger}
    colLook {mustBeInteger}
end


[in_rows,in_cols] = size(image);
if in_rows < rowLook || in_cols < colLook
    error("The image size must be bigger than look numbers.");
end

global gars_mtlk_gpu_enable
if isempty(gars_mtlk_gpu_enable)
    gars_mtlk_gpu_enable = true;
end

if gars_mtlk_gpu_enable
    try
        result = MTLK_native(image, in_rows, in_cols, rowLook, colLook);
        return;
    catch e
        warning(e.identifier, ...
            "An error occurred when calling library, fallback to matlab\n" + ...
            "        Error message: %s", e.message);
        gars_mtlk_gpu_enable = false;
    end
end

result = MTLK_matlab(image, in_rows, in_cols, rowLook, colLook);

end


%---------- native function caller ----------%

function result = MTLK_native(image, in_rows, in_cols, row_look, col_look)

out_rows = fix(in_rows / row_look);
out_cols = fix(in_cols / col_look);

if ~isreal(image)
    if class(image) == "double"
        [errno1,rpart] = clib.gars.Multilookd(real(image), row_look, ...
            col_look, out_rows, out_cols);
        [errno2,ipart] = clib.gars.Multilookd(imag(image), row_look, ...
            col_look, out_rows, out_cols);
    else
        [errno1,rpart] = clib.gars.Multilookf(real(image), row_look, ...
            col_look, out_rows, out_cols);
        [errno2,ipart] = clib.gars.Multilookf(imag(image), row_look, ...
            col_look, out_rows, out_cols);
    end

    if errno1 ~= 0 || errno2 ~= 0
        error("error number %d and %d is returned.", errno1, errno2);
    end

    result = rpart + 1i * ipart;
else
    if class(image) == "double"
        [errno,result] = clib.gars.Multilookd(image, row_look, col_look, ...
            out_rows, out_cols);
    else
        [errno,result] = clib.gars.Multilookf(image, row_look, col_look, ...
            out_rows, out_cols);
    end

    if errno ~= 0
        error("error number %d is returned.", errno);
    end
end

end


%---------- MATLAB function caller ----------%

function result = MTLK_matlab(image, in_rows, in_cols, row_look, col_look)

row_strides = repmat(row_look, 1, floor(in_rows / row_look));
row_left = rem(in_rows, row_look);

col_strides = repmat(col_look, 1, floor(in_cols / col_look));
col_left = rem(in_cols, col_look);

patched = mat2cell(image(1:end-row_left,1:end-col_left), ...
    row_strides, col_strides);

result = cellfun(@(x) mean(x(:)), patched);

end
