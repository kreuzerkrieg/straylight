cmake_minimum_required(VERSION 3.10)
project(StrayLight)

set(CMAKE_CXX_STANDARD 17)
set(SOURCES
        src/main.cpp
        src/config.h
        src/configStraylight.h
        src/forward_model_single_psf_dual_resolution.cpp
        src/forward_model_single_psf_dual_resolution.h
        src/harvey_brdf.cpp
        src/harvey_brdf.h
        src/harvey_psf.cpp
        src/harvey_psf.h
        src/particulate_contamination_harvey_psf.cpp
        src/particulate_contamination_harvey_psf.h
        src/tis_surface_scattering_harvey.cpp
        src/tis_surface_scattering_harvey.h
        src/utilities.cpp
        src/utilities.h)
add_compile_options(
        -fopenmp
)
add_definitions(-DUSE_DOUBLEPRECISION)
add_executable(StrayLight
        ${SOURCES}
        )

SET_TARGET_PROPERTIES(StrayLight
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
        )

TARGET_LINK_LIBRARIES(StrayLight
        fftw3
        -fopenmp
        rt
        )

#target_compile_definitions(StrayLight
#        PRIVATE USE_OPENMP=1
#        )

find_package(CUDA REQUIRED)

cuda_add_executable(StrayLightFull
        ${SOURCES}
        src/cudaConvolution.cu
        )

SET_TARGET_PROPERTIES(StrayLightFull
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
        )
target_compile_definitions(StrayLightFull
        PRIVATE USE_EIGEN=1
        PRIVATE USE_OPENMP=1
        PRIVATE USE_CUDA_GPU=1
        )
target_include_directories(StrayLightFull
        PRIVATE ${CUDA_INCLUDE_DIRS}
        )
TARGET_LINK_LIBRARIES(StrayLightFull
        fftw3
        ${CUDA_LIBRARIES}
        -fopenmp
        rt
        )

set(MODERNIZED_SOURCES
        src/Modernize/main.cpp
        src/Modernize/config.h
        src/Modernize/configStraylight.h
        src/Modernize/forward_model_single_psf_dual_resolution.cpp
        src/Modernize/forward_model_single_psf_dual_resolution.h
        src/Modernize/harvey_brdf.cpp
        src/Modernize/harvey_brdf.h
        src/Modernize/harvey_psf.cpp
        src/Modernize/harvey_psf.h
        src/Modernize/particulate_contamination_harvey_psf.cpp
        src/Modernize/particulate_contamination_harvey_psf.h
        src/Modernize/tis_surface_scattering_harvey.cpp
        src/Modernize/tis_surface_scattering_harvey.h
        src/Modernize/utilities.cpp
        src/Modernize/utilities.h
        src/Modernize/Matrices.cpp
        src/Modernize/Matrices.h
        )
add_definitions(-DUSE_DOUBLEPRECISION)
add_executable(ModernStrayLight
        ${MODERNIZED_SOURCES}
        )

SET_TARGET_PROPERTIES(ModernStrayLight
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
        )

TARGET_LINK_LIBRARIES(ModernStrayLight
        -fopenmp
        fftw3
        rt
        )
target_compile_options(ModernStrayLight
        #        PUBLIC -O3
        #        PUBLIC -ftree-vectorizer-verbose=2
        ##        PUBLIC -fopt-info
        #        PUBLIC -ftree-vectorize
        #        PUBLIC -fopt-info-vec-missed
        #        PUBLIC -mavx
        #        PUBLIC -floop-parallelize-all
        #        PUBLIC -ftree-parallelize-loops=4
        #        PUBLIC -ftree-slp-vectorize
        PUBLIC -Wall
        PUBLIC -Wextra
        #        PUBLIC -Weverything
        #        PUBLIC -Wno-c++98-compat
        #        -Wno-c++98-compat-pedantic
        PUBLIC -pedantic
        PUBLIC -O3
        PUBLIC -ftree-vectorize
        #        PUBLIC -ftree-vectorizer-verbose=2
        #        PUBLIC -floop-parallelize-all
        #        PUBLIC -ftree-parallelize-loops=4
        PUBLIC -march=native
        #        PUBLIC -fopt-info-vec-missed
        PUBLIC -msse
        PUBLIC -msse2
        PUBLIC -mavx
        PUBLIC -mavx2
        -fopenmp
        #        -Rpass=loop-vectorize
        #        -Rpass-missed=loop-vectorize
        #        -Rpass-analysis=loop-vectorize
        )
