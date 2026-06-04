% internal__Multilook_gpu - Internal GPU-accelerated implementation of multi-looking
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function result = internal__Multilook_gpu(image, row_look, col_look)

if ~ismember(class(image), ["double" "single" "int32"])
    error("The GPU version of 'Multilook' only supports 'double', " + ...
        "'single' and 'int32'.");
end

global IPSARTMexHost;
if isempty(IPSARTMexHost)
    IPSARTMexHost = mexhost();
end

if ~isreal(image)
    r_part = IPSARTMexHost.feval("internal__mex_bridge", "Multilook", ...
        real(image), int32(row_look), int32(col_look));
    i_part = IPSARTMexHost.feval("internal__mex_bridge", "Multilook", ...
        imag(image), int32(row_look), int32(col_look));
    result = r_part + 1i * i_part;
else
    result = IPSARTMexHost.feval("internal__mex_bridge", "Multilook", ...
        image, int32(row_look), int32(col_look));
end

end