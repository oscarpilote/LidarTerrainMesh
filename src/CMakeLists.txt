include_directories(
	${swiss_lidar_SOURCE_DIR}/include/
	${swiss_lidar_SOURCE_DIR}/extern/
	)

#[===[
add_library(miniply
	../extern/miniply/miniply.cpp
	)
#[===[

#[===[
add_library(meshoptimizer
	../extern/meshoptimizer/src/allocator.cpp
	../extern/meshoptimizer/src/vfetchoptimizer.cpp
	../extern/meshoptimizer/src/vcacheoptimizer.cpp
	../extern/meshoptimizer/src/indexcodec.cpp
	../extern/meshoptimizer/src/vertexcodec.cpp
	../extern/meshoptimizer/src/simplifier.cpp
	)
#[===[

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

add_executable(swiss_stitch
	mesh.cpp
	mesh_adjacency.cpp
	mesh_ply.cpp
	swiss_stitch.cpp
	)
target_link_libraries(swiss_stitch miniply)
target_compile_features(swiss_stitch PRIVATE cxx_std_17)
target_compile_options(swiss_stitch PRIVATE -Wall)

#[===[
add_executable(swiss_build_one
	mesh.cpp
	mesh_ply.cpp
	swiss_mgrid.cpp
	swiss_build_one.cpp
	)
target_link_libraries(swiss_build_one miniply meshoptimizer)
target_compile_features(swiss_build_one PRIVATE cxx_std_17)
target_compile_options(swiss_build_one PRIVATE -Wall)
]===]
#

add_executable(french_lidar
	chrono.cpp
	las_read.cpp
	las_source.cpp
	las_normal.cpp
	copc.cpp
	mesh.cpp
	mesh_adjacency.cpp
	mesh_ply.cpp
	mesh_inria.cpp
	mesh_remesh.cpp
	vertex_remap.cpp
	mesh_utils.cpp
	french_lidar.cpp
	)
target_compile_features(french_lidar PRIVATE cxx_std_17)
target_compile_options(french_lidar PRIVATE -Wall)
target_link_libraries(french_lidar miniply meshoptimizer lazperf ${swiss_lidar_SOURCE_DIR}/extern/mmg/libmmgs.a)

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

add_executable(ply_to_translated_obj
	mesh.cpp
	mesh_ply.cpp
	mesh_obj.cpp
	ply_to_translated_obj.cpp
	)
target_link_libraries(ply_to_translated_obj miniply)
target_compile_features(ply_to_translated_obj PRIVATE cxx_std_17)
target_compile_options(ply_to_translated_obj PRIVATE -Wall)

