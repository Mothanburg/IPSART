% CloudePottier - Cloude-Pottier H/α/A 分解
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function [H,alpha,A] = CloudePottier(M)

arguments
    M {mustBeA(M, ["PolM2" "PolM3"])}
end

try
    [H,alpha,A] = internal__CloudePottier_native(M);
    return;
catch e
    warning(e.identifier, ...
        "An error occurred when calling library, fallback to matlab. " + ...
        "Error message: \n%s", e.message);
    global IPSARTMexHost
    clear global IPSARTMexHost
    [H,alpha,A] = internal__CloudePottier_matlab(M);
end

end
