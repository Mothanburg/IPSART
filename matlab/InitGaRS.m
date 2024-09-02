function InitGaRS(options)

arguments
    options.imports(1,:) string = ["Utils", "Basic", "Data", "SAR", "PolSAR", "InSAR"]
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
% Check whether there is a manifest file or not.
manifest = fullfile(pkgpath, "manifest.txt");
if exist(manifest, "file")
    internals = readlines(manifest);
    for internal = internals'
        internal_path = fullfile(pkgpath, internal);
        if ~exist(internal_path, "dir")
            error("The internal package ""%s"" of ""%s"" doesn't exist, " + ...
                "the GaRS library may be broken.", internal, pkg);
        end
        addpath(internal_path);
    end
end

end

