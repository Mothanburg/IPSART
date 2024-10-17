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
                        os.vcp(libpath, installdir)
                        installed[libname] = true
                    end
                end
            end
        end

        -- install GaRS library file
        local tgt_file = target:targetfile()
        os.vcp(tgt_file, installdir)
        if target:is_plat("windows") then
            tgt_lib_file = string.gsub(tgt_file, "%.dll$", ".lib")
            os.vcp(tgt_lib_file, installdir)
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
