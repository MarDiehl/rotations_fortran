file(GLOB sources *.f90 *.c)

add_executable(test_rotations ${sources})
install (TARGETS test_rotations RUNTIME DESTINATION bin)
      
