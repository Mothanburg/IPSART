% RefinedLeeFilter - 精致Lee滤波
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
%
% Refined Lee filter with 7x7 window
function result = RefinedLeeFilting(M, lookNum)
arguments
    M {mustBeA(M, ["PolM2" "PolM3"])}
    lookNum {mustBeGreaterThanOrEqual(lookNum, 1)}
end

try
    result = internal__RefinedLeeFilting_gpu(M, lookNum);
    return;
catch e
    warning(e.identifier, ...
        "An error occurred when calling library, fallback to matlab. " + ...
        "Error message: \n%s", e.message);
    global IPSARTMexHost
    clear global IPSARTMexHost
    result = internal__RefinedLeeFilting_matlab(M, lookNum);
end

end
