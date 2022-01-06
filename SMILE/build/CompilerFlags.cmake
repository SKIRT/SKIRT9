# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# adjust C++ compiler flags for current target to our needs
# ------------------------------------------------------------------

set_property(TARGET ${TARGET} PROPERTY CXX_STANDARD 14)
set_property(TARGET ${TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(${TARGET} PRIVATE -Wall -W -pedantic)
    if (NO_EXTRA_WARNINGS)
        target_compile_options(${TARGET} PRIVATE -Wno-unused-parameter)
    else()
        target_compile_options(${TARGET} PRIVATE
            -Wdeprecated -Wextra-semi -Wold-style-cast -Wdouble-promotion
            -Wunused-exception-parameter -Wmissing-variable-declarations
            -Wconditional-uninitialized -Wswitch-enum -Wcovered-switch-default)
            # -Wconversion (ignore size_t to/from int conversions)
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(${TARGET} PRIVATE -Wall -W -pedantic)
    if (NO_EXTRA_WARNINGS)
        target_compile_options(${TARGET} PRIVATE -Wno-misleading-indentation -Wno-unused-parameter)
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(${TARGET} PRIVATE -fp-model precise -Wall)
    if (NO_EXTRA_WARNINGS)
        target_compile_options(${TARGET} PRIVATE -Wno-deprecated)
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(${TARGET} PRIVATE /wd4267 /wd4244)  # ignore size_t to/from int conversions
    if (NO_EXTRA_WARNINGS)
        target_compile_options(${TARGET} PRIVATE /wd4996)  # ignore unsafe C-style std functions
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
