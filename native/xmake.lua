add_rules("mode.debug", "mode.release")
set_languages("c++20", "c99")
set_warnings("all")

add_requires("eigen", "openblas", "openmp", "opencl")

target("GaRS")
    set_kind("shared")
    add_includedirs("include/", { public = true })
    add_files("src/*.cpp")
    add_packages("eigen", "openblas", "openmp", "opencl")
    add_defines("COMPILING_GARS")
    if is_mode("release") then
        add_defines("EIGEN_NO_DEBUG")
    end
    on_config("on_config")
    on_install("install")


target("test_runable")
    set_default(false)
    set_kind("binary")
    add_files("test/test_runable.cpp")
    add_tests("default")
    add_deps("GaRS")

