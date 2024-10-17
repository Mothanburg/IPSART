function InitGaRS(options)

arguments
    options.imports (1,:) string = [
        "Utils", ...
        "Basic", ...
        "Data", ...
        "SAR", ...
        "PolSAR", ...
        "InSAR", ...
        "Geometric"...
        ]
    options.debug = true
end

% Create the basic config of the GaRS library.
[root,~,~] = fileparts(mfilename("fullpath"));

% Import function packages.
try
    for pkg = options.imports
        import_package(root, pkg);
    end
catch e
    rehash("path");
    rethrow(e);
end

% Using out-of-process mode for safety
addpath(fullfile(root, "bin"));
if ~isMATLABReleaseOlderThan("R2023a") && options.debug
    clibConfiguration("gars", "ExecutionMode", "outofprocess");
end

% ------------- You can ONLY edit the part below ------------- %

% Set function alias
assignin("base", "LS1PH", @(x) HistStretch(x, "Linear Percent", 0, 99));
assignin("base", "LS2PH", @(x) HistStretch(x, "Linear Percent", 0, 98));
assignin("base", "LS5PH", @(x) HistStretch(x, "Linear Percent", 0, 95));
assignin("base", "LSOPT", @(x) HistStretch(x, "Optimized Linear"));

end


function import_package(root, pkg)

pkgpath = fullfile(root, pkg);
if ~exist(pkgpath, "dir")
    error("Package ""%s"" doesn't exist, the GaRS library may be broken.", pkg);
end

addpath(pkgpath);

% Add sub packages
pkg_files = dir(pkgpath);

for file = pkg_files'
    if file.isdir && file.name ~= "." && file.name ~= ".."
        addpath(fullfile(file.folder, file.name));
    end
end

end

