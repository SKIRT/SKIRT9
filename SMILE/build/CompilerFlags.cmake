# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# adjust C++ compiler flags for current target to our needs
# ------------------------------------------------------------------

set_property(TARGET ${TARGET} PROPERTY CXX_STANDARD 14)
set_property(TARGET ${TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|IntelLLVM")  # the Intel oneAPI compiler supports Clang options
    target_compile_options(${TARGET} PRIVATE -Wall -W -pedantic)
    if (NO_EXTRA_WARNINGS)
        target_compile_options(${TARGET} PRIVATE -Wno-unused-parameter -Wno-unused-function -Wno-sign-compare
            -Wno-deprecated-declarations -Wno-unused-variable -Wno-unused-but-set-variable
            -Wno-deprecated-copy-with-user-provided-copy)
    else()
        target_compile_options(${TARGET} PRIVATE
            -Wdeprecated -Wextra-semi -Wold-style-cast -Wdouble-promotion
            -Wunused-exception-parameter -Wmissing-variable-declarations
            -Wconditional-uninitialized -Wswitch-enum -Wcovered-switch-default)
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(${TARGET} PRIVATE -Wall -W -pedantic)
    if (NO_EXTRA_WARNINGS)
        target_compile_options(${TARGET} PRIVATE -Wno-misleading-indentation -Wno-unused-parameter
            -Wno-unused-function -Wno-unused-result -Wno-deprecated-copy -Wno-sign-compare -Wno-restrict
            -Wno-unused-variable -Wno-unused-but-set-variable -Wno-maybe-uninitialized -Wno-format)
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")  # this catches the deprecated Intel classic compiler
    target_compile_options(${TARGET} PRIVATE -fp-model precise -Wall)
    if (NO_EXTRA_WARNINGS)
        target_compile_options(${TARGET} PRIVATE -Wno-deprecated)
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(${TARGET} PRIVATE /wd4267 /wd4244)  # ignore size_t to/from int conversions
    if (NO_EXTRA_WARNINGS)
        target_compile_options(${TARGET} PRIVATE /wd2220 /wd4018 /wd4101 /wd4477 /wd4996)
    endif()
endif()

if (WARNINGS_AS_ERRORS)
    if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        target_compile_options(${TARGET} PRIVATE /WX)
    else()
        target_compile_options(${TARGET} PRIVATE -Werror)
    endif()
endif()

# ------------------------------------------------------------------
