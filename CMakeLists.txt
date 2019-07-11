PROJECT(CORK)

enable_language(C CXX)

cmake_minimum_required(VERSION 3.00)
set(CMAKE_CXX_STANDARD 11)

OPTION(USE_MINI_GMP "Use mini-gmp instead of gmp or mprr" ON)

IF(USE_MINI_GMP)
    add_definitions(-DUSE_MINI_GMP)

    include_directories(contrib/gmp-6.1.2/mini-gmp)
    
    SET(MINI_GMP_SOURCES    contrib/gmp-6.1.2/mini-gmp/mini-gmp.c
        contrib/gmp-6.1.2/mini-gmp/mini-gmp.h
        contrib/gmp-6.1.2/mini-gmp/mini-gmpxx.h
    )
ELSE(USE_MINI_GMP)
    SET(MINI_GMP_SRC )

    # find GMP library
    find_path(GMP_INCLUDE_DIR gmp.h gmpxx.h 
        PATHS /usr/include /usr/local/include /opt/local/include )

    find_library(GMP_LIBRARY NAMES gmp libgmp 
        PATHS /usr/lib /usr/local/lib /opt/local/lib)

    include_directories(${GMP_INCLUDE_DIR})
ENDIF(USE_MINI_GMP)

include_directories(src src/file_formats src/mesh src/util src/math src/rawmesh src/isct src/accel)

SET_SOURCE_FILES_PROPERTIES( src/isct/triangle.c PROPERTIES COMPILE_FLAGS "-Wno-format -DNO_TIMER -DREDUCED -DCDT_ONLY -DTRILIBRARY -DANSI_DECLARATORS" )


SET(CORK_SOURCES 
    src/file_formats/files.cpp
    src/file_formats/ifs.cpp
    src/file_formats/off.cpp
    src/util/log.cpp
    src/util/timer.cpp
    src/cork.cpp
    src/main.cpp
    src/isct/empty3d.cpp
    src/isct/triangle.c
    src/isct/quantization.cpp
)

SET(CORK_HEADERS
    src/file_formats/files.h
    src/mesh/mesh.h
    src/mesh/mesh.decl.h
    src/cork.h
    src/rawmesh/rawMesh.h
    src/util/iterPool.h
    src/util/shortVec.h
    src/util/unionFind.h
    src/util/memPool.h
    src/util/prelude.h
    src/math/bbox.h
    src/math/vec.h
    src/math/ray.h
    src/accel/aabvh.h
    src/isct/fixext4.h
    src/isct/quantization.h
    src/isct/ext4.h
    src/isct/gmpext4.h
    src/isct/triangle.h
    src/isct/empty3d.h
    src/isct/absext4.h
    src/isct/unsafeRayTriIsct.h
    src/isct/fixint.h
)

ADD_EXECUTABLE(cork 
    ${CORK_SOURCES} ${CORK_HEADERS} ${MINI_GMP_SOURCES}
)

#ADD_EXECUTABLE(test-cork 
#    src/test_main.cpp
#)

#TARGET_LINK_LIBRARIES(cork triangle)


ADD_EXECUTABLE(off2obj
    src/off2obj.cpp
)

#TARGET_LINK_LIBRARIES(msh_test msh_io rply)

IF(NOT USE_MINI_GMP)
    TARGET_LINK_LIBRARIES(cork ${GMP_LIBRARY})
ENDIF(NOT USE_MINI_GMP)
