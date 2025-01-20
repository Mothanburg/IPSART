% General four component decomposition with unitary matrix transformation
function [Ps,Pd,Pv,Ph] = G4U(T3)

arguments
    T3 PolT3
end

global MATSAR_G4U_ENABLE_CPU
if isempty(MATSAR_G4U_ENABLE_CPU)
    MATSAR_G4U_ENABLE_CPU = true;
end

if MATSAR_G4U_ENABLE_CPU
    try
        [Ps,Pd,Pv,Ph] = internal__G4U_native(T3);
        return;
    catch e
        warning(e.identifier, ...
            "An error occurred when calling library, fallback to matlab\n" + ...
            "        Error message: %s", e.message);
        MATSAR_G4U_ENABLE_CPU = false;
    end
end

[Ps,Pd,Pv,Ph] = internal__G4U_matlab(T3);

end
