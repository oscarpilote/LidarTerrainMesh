#/bin/bash
if [ $# -ge 1 ] && [ "$1" = "debug" ]
then
	echo "Building DEBUG version."
	echo "-----------------------"
	make --no-print-directory -C build/debug -j 8
else
	echo "Building RELEASE version."
	echo "-----------------------"
	make --no-print-directory -C build/release -j 8
fi
