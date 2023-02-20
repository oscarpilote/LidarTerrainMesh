echo "Building DEBUG version."
echo "-----------------------"
make --no-print-directory -C build/debug -j 4

echo
echo "Building RELEASE version."
echo "-----------------------"
make --no-print-directory -C build/release -j 4
