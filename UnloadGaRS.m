function UnloadGaRS

clear global;
if ~isMATLABReleaseOlderThan("R2023a")
    clibConfiguration("gars").unload()
end

rehash("path");

end