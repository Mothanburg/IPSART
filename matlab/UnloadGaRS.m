function UnloadGaRS

if libisloaded("GaRS")
    unloadlibrary("GaRS");
end

rehash("path");

clear global;

end