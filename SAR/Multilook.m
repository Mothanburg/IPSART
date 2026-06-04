% Multilook - SAR图像多视处理
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
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

try
    result = internal__Multilook_gpu(image, rowLook, colLook);
    return;
catch e
    warning(e.identifier, ...
        "An error occurred when calling library, fallback to matlab. " + ...
        "Error message: \n%s", e.message);
    global IPSARTMexHost
    clear global IPSARTMexHost
    result = internal__Multilook_matlab(image, rowLook, colLook);
end


end
