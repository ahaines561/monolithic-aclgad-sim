cmake_minimum_required(VERSION 3.10)
project(silvaco_parse)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(Garfield REQUIRED)
find_package(OpenMP REQUIRED)
# add_executable(test_silvaco test_silvaco.cpp)
# target_link_libraries(test_silvaco Garfield::Garfield)
add_executable(test_silvaco_Efield test_silvaco_Efield.cpp)
target_link_libraries(test_silvaco_Efield Garfield::Garfield)
add_executable(probe_mesh probe_mesh.cpp)
target_link_libraries(probe_mesh Garfield::Garfield)
add_executable(main main.cpp lgad_tools.cpp)
target_link_libraries(main Garfield::Garfield OpenMP::OpenMP_CXX)
add_executable(probe_field probe_field.cpp)
target_link_libraries(probe_field Garfield::Garfield)
# add_executable(plot_avalanche plot_avalanche.cpp)
# target_link_libraries(plot_avalanche Garfield::Garfield)
# add_executable(planar_movie planar_movie.cpp)
# target_link_libraries(planar_movie Garfield::Garfield)
# add_executable(attempt attempt.cpp)
# target_link_libraries(attempt Garfield::Garfield)
# add_executable(1000MIPs 1000MIPs.cpp lgad_tools.cpp)
# target_link_libraries(1000MIPs Garfield::Garfield OpenMP::OpenMP_CXX)
add_executable(feedback_scan feedback_scan.cpp lgad_tools.cpp)
target_link_libraries(feedback_scan Garfield::Garfield OpenMP::OpenMP_CXX)

# --- M2/M5 tools -------------------------------------------------------------
# sheet_psi_dump and sheet_prodcheck need only sheet_weighting.hh (header-only,
# no Garfield). garfield_transport needs MediumSilicon. See tools/README.md.
add_executable(sheet_psi_dump tools/sheet_psi_dump.cpp)
target_include_directories(sheet_psi_dump PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(sheet_prodcheck tools/sheet_prodcheck.cpp)
target_include_directories(sheet_prodcheck PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Closed-form check of the sheet solver against Riegler arXiv:1602.07949 Eq.138.
# The only EXTERNAL test of sheet_weighting.hh; everything else is self-consistency.
add_executable(riegler_benchmark tools/riegler_benchmark.cpp)
target_include_directories(riegler_benchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(garfield_transport tools/garfield_transport.cpp)
target_link_libraries(garfield_transport Garfield::Garfield)
# add_executable(rng_test rng_test.cpp lgad_tools.cpp)
# target_link_libraries(rng_test Garfield::Garfield OpenMP::OpenMP_CXX)