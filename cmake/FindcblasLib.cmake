# Find the cblas libraries as part of Accelerate.framework or as standalon framework
#
# The following are set after configuration is done:
#  CBLAS_FOUND
#  CBLAS_INCLUDE_DIR
#  CBLAS_LINKER_LIBS


if(NOT APPLE)
  return()
endif()

#include_directories(/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers)

#System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers/

#set(__veclib_include_suffix "Frameworks/vecLib.framework/Versions/Current/Headers")

#set(__cblas_include_suffix "Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers/")

find_path(CBLAS_INCLUDE_DIR cblas.h
          DOC "cblas include directory"
          PATHS /System/Library/Frameworks/Accelerate.framework/Versions/Current/${__veclib_include_suffix}
                /System/Library/${__veclib_include_suffix}
                /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Headers/
                /System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers
          NO_DEFAULT_PATH)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(vecLib DEFAULT_MSG vecLib_INCLUDE_DIR)

if(VECLIB_FOUND)
  if(vecLib_INCLUDE_DIR MATCHES "^/System/Library/Frameworks/vecLib.framework.*")
    set(vecLib_LINKER_LIBS -lcblas "-framework vecLib")
    message(STATUS "Found standalone vecLib.framework")
  else()
    set(vecLib_LINKER_LIBS -lcblas "-framework Accelerate")
    message(STATUS "Found vecLib as part of Accelerate.framework")
  endif()

  mark_as_advanced(vecLib_INCLUDE_DIR)
endif()