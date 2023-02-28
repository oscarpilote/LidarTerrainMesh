include_directories(
	${swiss_lidar_SOURCE_DIR}/include/
	${swiss_lidar_SOURCE_DIR}/extern/
	)

add_library(miniply
	../extern/miniply/miniply.cpp
	)

add_executable(swiss_lidar
	las_read.cpp
	las_source.cpp
	las_normal.cpp
	mesh.cpp
	mesh_ply.cpp
	mesh_inria.cpp
	vertex_remap.cpp
	mesh_utils.cpp
	swiss_lidar.cpp
	)

target_compile_features(swiss_lidar PRIVATE cxx_std_17)
target_compile_options(swiss_lidar PRIVATE -Wall)
target_link_libraries(swiss_lidar miniply)

add_executable(compute_normals
	las_read.cpp
	las_source.cpp
	las_normal.cpp
	mesh.cpp
	mesh_ply.cpp
	compute_normals.cpp
	)
target_link_libraries(compute_normals miniply)
target_compile_features(compute_normals PRIVATE cxx_std_17)
target_compile_options(compute_normals PRIVATE -Wall)
