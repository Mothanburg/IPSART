-- install GaRS native library
function main(target)

        if not is_mode("release", "releasedbg", "minsizerel") then
            local curmode = get_config("mode")
            raise("please change mode '%s' to a kind of release mode.", curmode)
        end

        local installdir = path.join("..", "matlab", "bin")
        installdir = path.absolute(installdir, os.projectdir())
        if not os.exists(installdir) then
            os.mkdir(installdir)
        end

        -- install GaRS.h
        local headerpath = path.join(os.projectdir(), "include", "GaRS.h")
        os.cp(headerpath, installdir)

        -- install all shared libs of depended packages
        local installed = {}
        for _, pkg in ipairs(target:orderpkgs()) do
            if pkg:enabled() then
                for _, libpath in ipairs(table.wrap(pkg:get("libfiles"))) do
                    if _is_shared_lib(target, libpath) then
                        local libname = path.filename(libpath)
                        if installed[libname] then
                            wprint("'%s' already exists, overwriting it.", libname)
                        end
                        os.cp(libpath, installdir)
                        installed[libname] = true
                    end
                end
            end
        end

        -- install GaRS.dll
        local tgt_file = target:targetfile()
        os.cp(tgt_file, installdir)

        -- if the platform is Windows, copy GaRS.lib for generating MATLAB interface
        if target:is_plat("windows") then
            os.cp((string.gsub(tgt_file, "%.dll", ".lib")), installdir)
        end
end

function _is_shared_lib(target, libpath)
    if target:is_plat("windows") then
        return libpath:endswith(".dll")
    else -- plat == "linux" or plat == "macosx" or any else
        return libpath:endswith(".so") or libpath:match(".+%.so%..+$") or libpath:endswith(".dylib")
    end
end