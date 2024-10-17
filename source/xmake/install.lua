import("core.project.config")

-- install GaRS native library
function main(target)

        local installdir = path.join(os.projectdir(), "..", "bin")
        if not os.exists(installdir) then
            os.mkdir(installdir)
        end

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

        -- install GaRS library file
        local tgt_file = target:targetfile()
        os.cp(tgt_file, installdir)
        if target:is_plat("windows") then
            tgt_lib_file = string.gsub(tgt_file, "%.dll$", ".lib")
            os.cp(tgt_lib_file, installdir)
        end

        -- for MSYS2, install the necessary libstdc++/libomp/libpthread/... runtime libraries
        if target:is_plat("msys") then
            local envs = os.getenvs()
            local msystem = envs["MSYSTEM"]
            local bin_dir = "/" .. string.lower(msystem) .. "/bin/"
            local runtimes_libs = {"libgcc_s_seh-1.dll", "libgomp-1.dll", "libwinpthread-1.dll", "libstdc++-6.dll"}
            for _, libname in pairs(runtimes_libs) do
                os.runv("cp", {bin_dir .. libname, installdir})
            end
            -- remove GaRS.lib to help MATLAB finding the correct library
            if os.exists(path.join(installdir, "GaRS.lib")) then
                os.rm(path.join(installdir, "GaRS.lib"))
            end
        end

        print("You may need to run the bootstrap script 'generate_matlab_interface.m' in the '../bootstrap'")

end

function _is_shared_lib(target, libpath)
    if target:is_plat("windows") then
        return libpath:endswith(".dll")
    else -- plat == "linux" or plat == "macosx" or any else
        return libpath:endswith(".so") or libpath:match(".+%.so%..+$") or libpath:endswith(".dylib")
    end
end
