% Cloude-Pottier H/a/A decomposition for dual-pol
function [H,alpha,A] = DualPolCloudePottier(M2)

arguments
    M2 PolM2
end

global MATSAR_DUALPOLCLOUDEPOTTIER_ENABLE_CPU
if isempty(MATSAR_DUALPOLCLOUDEPOTTIER_ENABLE_CPU)
    MATSAR_DUALPOLCLOUDEPOTTIER_ENABLE_CPU = true;
end

if MATSAR_DUALPOLCLOUDEPOTTIER_ENABLE_CPU
    try
        [H,alpha,A] = internal__DualPolCloudePottier_native(M2);
        return;
    catch e
        warning(e.identifier, ...
            "An error occurred when calling library, fallback to matlab\n" + ...
            "        Error message: %s", e.message);
        MATSAR_DUALPOLCLOUDEPOTTIER_ENABLE_CPU = false;
    end
end

[H,alpha,A] = internal__DualPolCloudePottier_matlab(M2);

end
