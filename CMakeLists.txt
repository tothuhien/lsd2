##################################################################
# LSD2 cmake build definition
# Added by Bui Quang Minh in April 2020
##################################################################

# Build Release mode by default with compiler optimizations turned on (e.g. -O3)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if (USE_LSD2)
    # compile LSD2 as a library
    add_library(lsd2
        src/confidence_interval.cpp
        src/confidence_interval.h
        src/date.h
        src/dating.cpp
        src/dating.h
        src/estimate_root.cpp
        src/estimate_root.h
        src/lsd.cpp
        src/lsd.h
        src/node.h
        src/options.cpp
        src/options.h
        src/outliers.cpp
        src/outliers.h
        src/pair.h
        src/part.h
        src/pr.h
        src/readData.cpp
        src/readData.h
        src/subtree.h
        src/utils.cpp
        src/utils.h
    )
else()
    # compile LSD2 as a stand-alone program
    add_executable(lsd2
        src/confidence_interval.cpp
        src/confidence_interval.h
        src/date.h
        src/dating.cpp
        src/dating.h
        src/estimate_root.cpp
        src/estimate_root.h
        src/lsd.cpp
        src/lsd.h
        src/node.h
        src/options.cpp
        src/options.h
        src/outliers.cpp
        src/outliers.h
        src/pair.h
        src/part.h
        src/pr.h
        src/readData.cpp
        src/readData.h
        src/subtree.h
        src/utils.cpp
        src/utils.h
    )
endif()

