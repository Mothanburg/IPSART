% Yamaguchi four-component decomposition (with negetive power elimination)
function [Ps,Pd,Pv,Ph] = Yamaguchi(T3)

arguments
    T3 PolT3
end

global MATSAR_YAMAGUCHI_ENABLE_CPU
if isempty(MATSAR_YAMAGUCHI_ENABLE_CPU)
    MATSAR_YAMAGUCHI_ENABLE_CPU = true;
end

if false%MATSAR_YAMAGUCHI_ENABLE_CPU
    try
        [Ps,Pd,Pv,Ph] = internal__Yamaguchi_native(T3);
        return;
    catch e
        warning(e.identifier, ...
            "An error occurred when calling library, fallback to matlab\n" + ...
            "        Error message: %s", e.message);
        MATSAR_YAMAGUCHI_ENABLE_CPU = false;
    end
end

[Ps,Pd,Pv,Ph] = internal__Yamaguchi_matlab(T3);

end
