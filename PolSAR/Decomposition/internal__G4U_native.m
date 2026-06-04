% internal__G4U_native - Internal MEX bridge for G4U decomposition
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function [Ps,Pd,Pv,Ph] = internal__G4U_native(T3)

global IPSARTMexHost;
if isempty(IPSARTMexHost)
    IPSARTMexHost = mexhost();
end

[Ps,Pd,Pv,Ph] = IPSARTMexHost.feval("internal__mex_bridge", "G4U", T3.m11, ...
    T3.m22, T3.m33, T3.m12_r, T3.m13_r, T3.m23_r, T3.m12_i, T3.m13_i, T3.m23_i);

end