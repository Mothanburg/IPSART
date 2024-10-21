import("core.project.config")

-- install GaRS native library
function main(target)

        local installdir = "../bin/"
        if not os.exists(installdir) then
            os.mkdir(installdir)
        else
            for _, existed in ipairs(os.files(installdir .. "*|*garsInterface*")) do
                os.rm(installdir .. existed)
            end
        end

        -- Install GaRS library file itself
        local tgt_file = target:targetfile()
        os.cp(tgt_file, installdir)

        -- Install the import library for MATLAB bootstrapping if compiled by MSVC
        if target:is_plat("windows") then
            local bootstrapdir = "../bootstrap/"
            local tgt_lib_file = string.gsub(tgt_file, "%.dll$", ".lib")
            os.cp(tgt_lib_file, bootstrapdir)
        end

        -- Install all dependency libraries for common platforms
        -- TODO: fix symbolic links on linux systems
        for _, pkg in ipairs(target:orderpkgs()) do
            if pkg:enabled() then
                for _, libpath in ipairs(table.wrap(pkg:get("libfiles"))) do
                    if _is_shared_lib(target, libpath) then
                        os.cp(libpath, installdir)
                    end
                end
            end
        end

        -- We must install all dependency libraries manually (stdc++/omp/pthread/...) for MSYS2
        --     because they are all 'system' libraries
        if target:is_plat("msys") then
            local msystem = string.lower(os.getenvs()["MSYSTEM"])
            local prefix = "/" .. msystem .. "/bin/"
            local dependencies, _ = os.iorun("ldd " .. target:targetfile())
            -- The result of 'ldd' include windows system dll, we don't need to copy them.
            for dll in string.gmatch(dependencies, "([^%s]+%.[dD][lL][lL]) =>") do
                -- All dependencies should be placed in "/${MSYSTEM}/bin/"
                local dll_file = prefix .. dll
                try {
                    function ()
                        os.runv("cp", {dll_file, installdir})
                    end
                }
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
