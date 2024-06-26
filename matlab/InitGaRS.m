function InitGaRS(options)

arguments
    options.imports(1,:) string = ["Basic", "Data", "Polarimetry"]
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
use_native = true;
if ispc
    nativelib = pjoin(root, "clib", "GaRS.dll");
    if ~exist(nativelib, "file")
        use_native = false;
    end
elseif isunix
    nativelib = pjoin(root, "clib", "GaRS.so");
    if ~exist(nativelib, "file")
        use_native = false;
    end
elseif ismac
    nativelib = pjoin(root, "clib", "GaRS.dylib");
    if ~exist(nativelib, "file")
        nativelib = pjoin(root, "clib", "GaRS.so");
        if ~exist(nativelib, "file")
            use_native = false;
        end
    end
else
    use_native = false;
end

% Load the native library.
if use_native
    if libisloaded("GaRS")
        warning("The status of the GaRS library may be unhealthy. " + ...
            "Perhaps you should run ""UnloadGaRS""");
        unloadlibrary("GaRS");
    end
    loadlibrary(nativelib);

    GARS_CONFIG.NativeFunctions = libfunctions("GaRS");

    % todo: check the platform of the native library.

    GARS_CONFIG.PlatformLevel = 1;
else
    GARS_CONFIG.PlatformLevel = 0;
end


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

pkgpath = pjoin(root, pkg);
if ~exist(pkgpath, "dir")
    error("Package ""%s"" doesn't exist, the GaRS library may be broken.", pkg);
end

addpath(pkgpath);
% Check whether there is a manifest file or not.
manifest = pjoin(pkgpath, "manifest.txt");
if exist(manifest, "file")
    internals = readlines(manifest);
    for internal = internals'
        internal_path = pjoin(pkgpath, internal);
        if ~exist(internal_path, "dir")
            error("The internal package ""%s"" of ""%s"" doesn't exist, " + ...
                "the GaRS library may be broken.", internal, pkg);
        end
        addpath(internal_path);
    end
end

end

function path = pjoin(varargin)
len = 2 * nargin - 1;
full_args = cell(1, len);
full_args(1:2:len) = varargin;
full_args(2:2:len) = {filesep};
path = strcat(full_args{:});
end

function loginfo(msg)
if evalin("caller", "options.verbose")
    disp(msg);
end
end
