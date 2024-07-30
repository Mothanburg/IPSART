function result = Multilook(image, rowLook, colLook)

arguments
    image (:,:)
    rowLook {mustBeInteger}
    colLook {mustBeInteger}
end


[in_rows,in_cols] = size(image);
if in_rows < rowLook || in_cols < colLook
    error("The image size must be bigger than look numbers.")
end

% We drop the last rows and cols whose length is less than look number.

global GARS_CONFIG

if GARS_CONFIG.CAPABILITY >= 1 && class(image) == "single" && numel(image) > 400000

    try
        if ~isreal(image)
            [errno,r_part] = clib.gars.Multilook(real(image), rowLook, colLook, ...
                fix(in_rows / rowLook), fix(in_cols / colLook));
            if errno ~= 0
                error("error number %d is returned.", errno);
            end

            [errno,i_part] = clib.gars.Multilook(imag(image), rowLook, colLook, ...
                fix(in_rows / rowLook), fix(in_cols / colLook));
            if errno ~= 0
                error("error number %d is returned.", errno);
            end
            
            result = r_part + 1i * i_part;
        else
            [errno,result] = clib.gars.Multilook(image, rowLook, colLook, ...
                fix(in_rows / rowLook), fix(in_cols / colLook));
            if errno ~= 0
                error("error number %d is returned.", errno);
            end
        end

        return;
    catch e
        warning(e.identifier, "An error occurred when calling library, fallback to matlab.\n" + ...
            "        Error message: %s", e.message);
    end
end

result = multlk_matlab(image, in_rows, in_cols, rowLook, colLook);

end


function result = multlk_matlab(image, in_rows, in_cols, rowLook, colLook)

row_strides = repmat(rowLook, 1, floor(in_rows / rowLook));
row_left = rem(in_rows, rowLook);

col_strides = repmat(colLook, 1, floor(in_cols / colLook));
col_left = rem(in_cols, colLook);

patched = mat2cell(image(1:end-row_left,1:end-col_left), row_strides, col_strides);

result = cellfun(@(x) mean(x(:)), patched);

end
