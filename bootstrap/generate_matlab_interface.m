%% Run this part to generate intermediate files

header = "../include/GaRS.h";

if ispc
    lib = "../bin/GaRS.lib";
else
    if exist("../bin/GaRS.dylib", "file")
        lib = "../bin/GaRS.dylib";
    else
        lib = "../bin/GaRS.so";
    end
end

clibgen.generateLibraryDefinition(...
    header, ...
    Libraries=lib, ...
    PackageName="gars", ...
    CLinkage=true, ...
    OutputFolder=".", ...
    Verbose=true, ...
    OverwriteExistingDefinitionFiles=true ...
);


%% You may check the 'definegars.m' before runing this part;

% Generate matlab interface library
build(library_definition());

copyfile("gars/garsInterface.dll", "../bin");

rmdir("gars", "s");
delete("definegars.m", "garsData.xml");