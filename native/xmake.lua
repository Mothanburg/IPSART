add_rules("mode.debug", "mode.release")

set_runtimes("MD")

set_warnings("all")

add_requires("eigen", "openblas", "openmp", "opencl")

target("GaRS")
    set_kind("shared")

    add_packages("eigen", "openblas", "openmp", "opencl")

    set_languages("c++20")
    add_includedirs("include/", { public = true })
    set_pcxxheader("include/pch.h")
    add_files("src/*.cpp")

    add_cxxflags("cl::/bigobj", "cl::/D_SILENCE_STDEXT_ARR_ITERS_DEPRECATION_WARNING")

    add_defines("COMPILING_GARS")
    if is_mode("release") then
        add_defines("EIGEN_NO_DEBUG")
    end

    set_installdir(path.join(os.projectdir(), "../matlab/bin"))
    on_config("xmake/config")
    on_install("xmake/install")

