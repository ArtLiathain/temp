# Docs
if ((${BUILD_DOC} STREQUAL "ONLY") OR (${BUILD_DOC} STREQUAL "YES"))
    message(STATUS "Requested to build the documentation")
    find_package(Doxygen 1.8 REQUIRED)
    find_package(LATEX REQUIRED COMPONENTS DVIPS)
elseif (${BUILD_DOC} STREQUAL "TRY")
    find_package(Doxygen 1.8)
    find_package(LATEX COMPONENTS DVIPS)
endif ()

if(DOXYGEN_FOUND AND LATEX_FOUND)
    # This could be updated to use doxygen_add_docs() rather than configuring existing Doxyfiles
    file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/doc/doc/html")
    configure_file("${PROJECT_SOURCE_DIR}/Doxygen.in"
                 "${PROJECT_BINARY_DIR}/Doxyfile" @ONLY)
    configure_file("${PROJECT_SOURCE_DIR}/DoxygenPyM.in"
                 "${PROJECT_BINARY_DIR}/DoxyfilePyM" @ONLY)
    add_custom_target(doc ALL Doxygen::doxygen
                        ${PROJECT_BINARY_DIR}/DoxyfilePyM
                        DEPENDS docrunner
                        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/doc/doc/html"
                        COMMENT "Generating documentation for developers"
                        VERBATIM)
    add_custom_target(docrunner ALL Doxygen::doxygen
                        ${PROJECT_BINARY_DIR}/Doxyfile
                        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/doc"
                        COMMENT "Generating documentation for developers"
                        VERBATIM)
    install(DIRECTORY "${PROJECT_BINARY_DIR}/doc/doc/html" DESTINATION
                share/ffea/doc OPTIONAL)
endif()
