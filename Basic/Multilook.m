function result = Multilook(image, rowLook, colLook)

arguments
    image (:,:)
    rowLook {mustBeInteger}
    colLook {mustBeInteger}
end

if rowLook == 1 && colLook == 1
    result = image;
    return;
end

[in_rows,in_cols] = size(image);
if in_rows < rowLook || in_cols < colLook
    error("The image size must be bigger than look numbers.");
end

global MATSAR_MULTILOOK_ENABLE_GPU
if isempty(MATSAR_MULTILOOK_ENABLE_GPU)
    MATSAR_MULTILOOK_ENABLE_GPU = true;
end

if MATSAR_MULTILOOK_ENABLE_GPU
    try
        result = internal__Multilook_gpu( ...
            image, in_rows, in_cols, rowLook, colLook);
        return;
    catch e
        warning(e.identifier, ...
            "An error occurred when calling library, fallback to matlab\n" + ...
            "        Error message: %s", e.message);
        MATSAR_MULTILOOK_ENABLE_GPU = false;
    end
end

result = internal__Multilook_matlab(image, in_rows, in_cols, rowLook, colLook);

end
