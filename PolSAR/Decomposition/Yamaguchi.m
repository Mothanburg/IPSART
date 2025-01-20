% Yamaguchi four-component decomposition
function [Ps,Pd,Pv,Ph] = Yamaguchi(C3)

arguments
    C3 PolC3
end

global MATSAR_YAMAGUCHI_ENABLE_CPU
if isempty(MATSAR_YAMAGUCHI_ENABLE_CPU)
    MATSAR_YAMAGUCHI_ENABLE_CPU = true;
end

if MATSAR_YAMAGUCHI_ENABLE_CPU
    try
        [Ps,Pd,Pv,Ph] = internal__Yamaguchi_native(C3);
        return;
    catch e
        warning(e.identifier, ...
            "An error occurred when calling library, fallback to matlab\n" + ...
            "        Error message: %s", e.message);
        MATSAR_YAMAGUCHI_ENABLE_CPU = false;
    end
end

[Ps,Pd,Pv,Ph] = internal__Yamaguchi_matlab(C3);

end
