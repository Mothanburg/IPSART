% internal__RefinedLeeFilter_gpu - Internal GPU-accelerated implementation of Refined Lee filter
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function result = internal__RefinedLeeFilting_gpu(M, look_num)

global IPSARTMexHost;
if isempty(IPSARTMexHost)
    IPSARTMexHost = mexhost();
end

if isa(M, "PolM2")
    [m11,m22,m12r,m12i] = IPSARTMexHost.feval("internal__mex_bridge", ...
        "RefinedLeeFilter", M.m11, M.m22, M.m12_r, M.m12_i, look_num);
    if isa(M, "PolC2")
        result = PolC2(m11, m22, m12r, m12i, M.PolType);
    elseif isa(M, "PolT2")
        result = PolT2(m11, m22, m12r, m12i);
    end
else
    [m11,m22,m33,m12r,m13r,m23r,m12i,m13i,m23i] = IPSARTMexHost.feval( ...
        "internal__mex_bridge", "RefinedLeeFilter", M.m11, M.m22, M.m33, ...
        M.m12_r, M.m13_r, M.m23_r, M.m12_i, M.m13_i, M.m23_i, look_num);
    if isa(M, "PolC3")
        result = PolC3(m11, m22, m33, m12r, m13r, m23r, m12i, m13i, m23i);
    elseif isa(M, "PolT3")
        result = PolT3(m11, m22, m33, m12r, m13r, m23r, m12i, m13i, m23i);
    end
end

end
