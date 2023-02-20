echo "Installing DEBUG build system into build/debug."
echo "-----------------------------------------------"
mkdir -p build/debug
cmake -B build/debug -DCMAKE_BUILD_TYPE=Debug

echo
echo "Installing RELEASE build system into build/debug."
echo "-----------------------------------------------"
mkdir -p build/release
cmake -B build/release -DCMAKE_BUILD_TYPE=Release

echo
./make.sh
