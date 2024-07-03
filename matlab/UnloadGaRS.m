function UnloadGaRS

rehash("path");

global GARS_CONFIG
if ~isMATLABReleaseOlderThan("R2023a")
    GARS_CONFIG.CLIB.unload();
end

clear GARS_CONFIG;

end