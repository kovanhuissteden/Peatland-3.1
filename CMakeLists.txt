project(peatland30)

add_executable(peatland30 main.cpp)

add_executable(peatland30 main.cpp matrix.cpp general.cpp readparams.cpp initialize.cpp paramcheck.cpp bookkeep.cpp heat.cpp water.cpp methane.cpp SOMproduction.cpp SOMdecomposition.cpp)

install(TARGETS peatland30 RUNTIME DESTINATION bin)
