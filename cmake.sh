#/bin/bash
if [ $# -ge 1 ] && [ "$1" = "debug" ]
then
	echo "Installing DEBUG build system into build/debug."
	echo "-----------------------------------------------"
	mkdir -p build/debug
	cmake -B build/debug -DCMAKE_BUILD_TYPE=Debug
	./make.sh $1
else
	echo "Installing RELEASE build system into build/debug."
	echo "-----------------------------------------------"
	mkdir -p build/release
	cmake -B build/release -DCMAKE_BUILD_TYPE=Release
	./make.sh
fi
