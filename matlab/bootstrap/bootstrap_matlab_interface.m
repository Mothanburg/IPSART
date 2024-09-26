%% You must run this part first.

header = "../../native/include/GaRS.h";

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
    Verbose=true ...
);

% Please Edit the define script by the instructions in the file.
edit(fullfile("build", "definegars.m"));

%% Then you should run this part.
addpath("H:\Code\GaRS\native\build");
build(definegars());
copyfile(fullfile("H:\Code\GaRS\native\build", "gars", "garsInterface.dll"), "H:\Code\GaRS\matlab\bin");
rehash("path");
