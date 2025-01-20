% Cloude-Pottier H/a/A decomposition
function [H,alpha,A] = CloudePottier(M3)

arguments
    M3 PolM3
end

if ~isa(M3, "PolT3")
    warning("The input data is not the coherency matrix T, " + ...
        "the returned alpha may not make sense");
end

global MATSAR_CLOUDEPOTTIER_ENABLE_CPU
if isempty(MATSAR_CLOUDEPOTTIER_ENABLE_CPU)
    MATSAR_CLOUDEPOTTIER_ENABLE_CPU = true;
end

if MATSAR_CLOUDEPOTTIER_ENABLE_CPU
    try
        [H,alpha,A] = internal__CloudePottier_native(M3);
        return;
    catch e
        warning(e.identifier, ...
            "An error occurred when calling library, fallback to matlab\n" + ...
            "        Error message: %s", e.message);
        MATSAR_CLOUDEPOTTIER_ENABLE_CPU = false;
    end
end

[H,alpha,A] = internal__CloudePottier_matlab(M3);

end
