echo "Cloning PS3D"
cd $SRC_DIR
git clone https://github.com/matt-frey/ps3d.git
cd ps3d
git checkout main-cheby-drew
git log | head

echo "Configuring PS3D"

./bootstrap
mkdir -p build
cd build
$SRC_DIR/ps3d/configure \
    --enable-verbose    \
    --enable-buoyancy   \
    --prefix=${PREFIX}/ps3d


echo "Compiling PS3D"
make

echo "Instaling PS3D"
make install
