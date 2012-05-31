SET (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
INCLUDE_DIRECTORIES (/usr/local/include/eigen3)
ADD_EXECUTABLE (qsAnnni main.cpp readinput.h readinput.cpp diag.h diag.cpp mps.h mps.cpp tsrproj.h tsrproj.cpp vmps.h vmps.cpp)
