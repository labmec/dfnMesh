# Code testing. Please run these before Pushing or PR
enable_testing()

# Simple
add_test(NAME Minimal COMMAND dfnGoBack "examples/minimal.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
# add_test(NAME Minimal2 COMMAND dfnGoBack "examples/minimal3Da.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME Tetrahedra_0 COMMAND dfnGoBack "examples/tetrahedra0.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Snapping and incorporating faces of the mesh
add_test(NAME Incorporate_1 COMMAND dfnGoBack "examples/incorporate_1.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME Incorporate_2 COMMAND dfnGoBack "examples/incorporate_2.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME Incorporate_3 COMMAND dfnGoBack "examples/incorporate_3.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME Incorporate_4 COMMAND dfnGoBack "examples/incorporate_4.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME Tetrahedra_1 COMMAND dfnGoBack "examples/tetrahedra1.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME Tetrahedra_2 COMMAND dfnGoBack "examples/tetrahedra2.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Limit recovery
add_test(NAME FracLimits1 COMMAND dfnGoBack "examples/fraclimits1.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME FracLimits2 COMMAND dfnGoBack "examples/fraclimits2.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME FracLimits3 COMMAND dfnGoBack "examples/fraclimits3.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME FracLimits4 COMMAND dfnGoBack "examples/fraclimits4.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME FracLimits5 COMMAND dfnGoBack "examples/fraclimits5.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME Octagon COMMAND dfnGoBack "examples/exampleOctagon.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Snap Induced Overlap
add_test(NAME SnapOverlap1 COMMAND dfnGoBack "examples/bug_snap_overlap1.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME SnapOverlap2 COMMAND dfnGoBack "examples/bug_snap_overlap2.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Splinter formation AKA fracture self-overlap
add_test(NAME SplinterFormation COMMAND dfnGoBack "examples/splinter1.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Frac 3D benchmarks
add_test(NAME Fl_Benchmark_3 COMMAND dfnGoBack "examples/flemisch_benchmark/fl_case3.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    # I ran into an infinite loop when logging DFNFaces. I fixed it, but here's a timeout last resource just to be safe. Ctest's default timeout is 1500 which is too damn long
    set_tests_properties(Fl_Benchmark_3 PROPERTIES TIMEOUT 20)
add_test(NAME Fl_Benchmark_2 COMMAND dfnGoBack "examples/flemisch_benchmark/fl_case2.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME Fl_Benchmark_1 COMMAND dfnGoBack "examples/flemisch_benchmark/fl_case1.json" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
# Benchmark4 is way too big to be a unit test