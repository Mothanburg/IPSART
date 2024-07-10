function InitGaRS(options)

arguments
    options.imports(1,:) string = ["Utils", "Basic", "Data", "Polarimetry"]
    options.verbose = false
end

% Create the basic config of the GaRS library.
global GARS_CONFIG;
GARS_CONFIG = struct();

[root,~,~] = fileparts(mfilename("fullpath"));

GARS_CONFIG.GaRSRoot = root;

% Import function packages.
try
    for pkg = options.imports
        import_package(root, pkg);
    end
catch e
    rehash("path");
    rethrow(e);
end

GARS_CONFIG.ImportedPackages = options.imports;

% Check the existence of the native library.
libpath = fullfile(root, "bin");
if exist(libpath, "dir")
    addpath(libpath);
    % For safety, we use out-of-process execution mode (if we can)
    if ~isMATLABReleaseOlderThan("R2023a")
        GARS_CONFIG.CLIB = clibConfiguration("gars", "ExecutionMode", "outofprocess");
    end
    % Test gpu capability
    try
        errno = clib.gars.GaRSTestOpenCL();
        if errno == 0
            loginfo("'%d' returned when testing GPU capability.");
            GARS_CONFIG.CAPABILITY = 2; % The platform can use GPU
        else
            GARS_CONFIG.CAPABILITY = 1; % The platform can only use CPU
        end
    catch
        GARS_CONFIG.CAPABILITY = 0;     % No native library capability
        rmpath(libpath);
    end
else
    GARS_CONFIG.CAPABILITY = 0;
end

% Start parallel pool
% GARS_CONFIG.POOL = parpool();


GARS_CONFIG.CAUTION= "This struct is crucial for the GaRS library, " + ...
    "DO NOT modify or clear it.";

% ------------- You can ONLY edit the part below ------------- %

% Set function alias.
assignin("base", "LS1PH", @(x) HistStretch(x, "Linear Percent", 0, 99));
assignin("base", "LS2PH", @(x) HistStretch(x, "Linear Percent", 0, 98));
assignin("base", "LS5PH", @(x) HistStretch(x, "Linear Percent", 0, 95));
assignin("base", "LSOPT", @(x) HistStretch(x, "Optimized Linear"));

% ------------- You can ONLY edit the part above ------------- %

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

function loginfo(msg, varargin)
if evalin("caller", "options.verbose")
    fmt = sprintf(msg, varargin{:});
    disp(fmt);
end
end
