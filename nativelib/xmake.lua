add_rules("mode.debug", "mode.release")
set_languages("c++20", "c99")
set_warnings("all")

add_requires("eigen", "openblas", "openmp")

target("cGaRS")
    set_kind("shared")
    add_includedirs("include/", { public = true })
    add_files("src/*.cpp")
    add_packages("eigen", "openblas", "openmp")
    add_defines("COMPILING_GARS")


target("test_cloudepottier")
    set_kind("binary")
    add_files("test/test_cloudepottier.cpp")
    add_tests("default")
    add_deps("cGaRS")

