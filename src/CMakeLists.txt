include_directories(
	${swiss_lidar_SOURCE_DIR}/include/
	${swiss_lidar_SOURCE_DIR}/extern/
	)

add_library(miniply
	../extern/miniply/miniply.cpp
	)

add_library(meshoptimizer
	../extern/meshoptimizer/src/allocator.cpp
	../extern/meshoptimizer/src/vfetchoptimizer.cpp
	../extern/meshoptimizer/src/vcacheoptimizer.cpp
	../extern/meshoptimizer/src/indexcodec.cpp
	../extern/meshoptimizer/src/vertexcodec.cpp
	../extern/meshoptimizer/src/simplifier.cpp
	)

add_executable(swiss_lidar
	chrono.cpp
	las_read.cpp
	las_source.cpp
	las_normal.cpp
	mesh.cpp
	mesh_adjacency.cpp
	mesh_ply.cpp
	mesh_inria.cpp
	mesh_remesh.cpp
	vertex_remap.cpp
	mesh_utils.cpp
	swiss_raster.cpp
	swiss_lidar.cpp
	)

target_compile_features(swiss_lidar PRIVATE cxx_std_17)
target_compile_options(swiss_lidar PRIVATE -Wall)
target_link_libraries(swiss_lidar miniply meshoptimizer gdal ${swiss_lidar_SOURCE_DIR}/extern/mmg/libmmgs.a)

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
