% internal__CloudePottier_native - Internal MEX bridge for Cloude-Pottier decomposition
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function [H,alpha,A] = internal__CloudePottier_native(M)

global IPSARTMexHost;
if isempty(IPSARTMexHost)
    IPSARTMexHost = mexhost();
end

if isa(M, "PolM2")
    [H,alpha,A] = IPSARTMexHost.feval("internal__mex_bridge", "CloudePottier", ...
        M.m11, M.m22, M.m12_r, M.m12_i);
else
    [H,alpha,A] = IPSARTMexHost.feval("internal__mex_bridge", "CloudePottier", ...
        M.m11, M.m22, M.m33, M.m12_r, M.m13_r, M.m23_r, M.m12_i, M.m13_i, ...
        M.m23_i);
end

end