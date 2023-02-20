include_directories(
	${swiss_lidar_SOURCE_DIR}/include/
	)

add_executable(swiss_lidar
	las_read.cpp
	las_source.cpp
	las_normal.cpp
	swiss_lidar.cpp
	)

target_compile_features(swiss_lidar PRIVATE cxx_std_17)
target_compile_options(swiss_lidar PRIVATE -Wall -Wextra)


