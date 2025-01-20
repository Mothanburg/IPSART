% Refined Lee filter with 7x7 window
function result = RefinedLeeFilting(M, lookNum)

arguments
    M PolMat
    lookNum {mustBeGreaterThanOrEqual(lookNum, 1)}
end

global MATSAR_REFINEDLEEFILTING_ENABLE_GPU;
if isempty(MATSAR_REFINEDLEEFILTING_ENABLE_GPU)
    MATSAR_REFINEDLEEFILTING_ENABLE_GPU = true;
end

if MATSAR_REFINEDLEEFILTING_ENABLE_GPU
    try
        result = internal__RefinedLeeFilting_gpu(M, lookNum);
        return;
    catch e
        warning(e.identifier, ...
            "An error occurred when calling library, fallback to matlab\n" + ...
            "        Error message: %s", e.message);
        MATSAR_REFINEDLEEFILTING_ENABLE_GPU = false;
    end
end

result = internal__RefinedLeeFilting_matlab(M, lookNum);

end
