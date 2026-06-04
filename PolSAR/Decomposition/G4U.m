% G4U - General four-component decomposition with unitary matrix transformation
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function [Ps,Pd,Pv,Ph] = G4U(T3)

arguments
    T3 PolT3
end

try
    [Ps,Pd,Pv,Ph] = internal__G4U_native(T3);
    return;
catch e
    warning(e.identifier, ...
        "An error occurred when calling library, fallback to matlab. " + ...
        "Error message: \n%s", e.message);
    global IPSARTMexHost
    clear global IPSARTMexHost
    [Ps,Pd,Pv,Ph] = internal__G4U_matlab(T3);
end



end
