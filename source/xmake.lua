add_rules("mode.debug", "mode.release")
set_defaultarchs("x64")
set_runtimes("MD")
set_warnings("all")

-- Dependencies of GaRS
dependencies = {
    "opencl",
    "openmp"
}
-- We only use 'system' packages (only installed by Pacman) in MSYS2 platform
if is_plat("msys") then
    add_requireconfs("*", { system = true })
    -- Workaround: Package 'Eigen' in MSYS2 is called 'eigen3' but not 'eigen', 
    --             so we must specifiy it explicitly
    table.insert(dependencies, "eigen3")
else
    -- Package 'Eigen' is called 'eigen' in xmake-repo
    table.insert(dependencies, "eigen")
end

add_requires(table.unpack(dependencies))

target("GaRS")
    set_kind("shared")
    set_prefixname("")
    add_packages(table.unpack(dependencies))

    set_languages("c++20")
    add_includedirs("../include", { public = true })
    add_includedirs("inc")
    set_pcxxheader("inc/pch.h")
    add_files("src/*.cpp")

    -- Enable necessary flags for MSVC
    add_cxxflags(
        "cl::/bigobj",
        "cl::/D_SILENCE_STDEXT_ARR_ITERS_DEPRECATION_WARNING"
        )

    add_defines("COMPILING_GARS")
    if is_mode("release") then
        add_defines("EIGEN_NO_DEBUG")
    end

    on_config("xmake/config")
    on_install("xmake/install")


target("testing")
    set_default(false)
    set_kind("binary")
    add_deps("GaRS")

    set_languages("c++20")

    add_files("test/main.cpp")

    on_install(function (target) 
        wprint("This target cannot be installed.")
    end)
