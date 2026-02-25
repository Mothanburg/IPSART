% Yamaguchi four-component decomposition (with negetive power elimination)
function [Ps,Pd,Pv,Ph] = Yamaguchi(T3)

arguments
    T3 PolT3
end

try
    [Ps,Pd,Pv,Ph] = internal__Yamaguchi_native(T3);
    return;
catch e
    warning(e.identifier, ...
        "An error occurred when calling library, fallback to matlab. " + ...
        "Error message: \n%s", e.message);
    global IPSARTMexHost
    clear global IPSARTMexHost
    [Ps,Pd,Pv,Ph] = internal__Yamaguchi_matlab(T3);
end

end
