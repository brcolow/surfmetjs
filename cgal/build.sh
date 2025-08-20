#! /bin/sh

# sudo apt-get install libcgal-dev libcgal-qt5-dev libgmp-dev libmpfr-dev

# Compile (typical):
# g++ -std=c++17 -O2 cgal.cpp -o cgal $(pkg-config --cflags --libs cgal) -lgmp -lmpfr

# Or without pkg-config (may vary by distro):
cgal="CGAL-6.1-beta1.tar.xz"
if ! [ -f "$cgal" ]; then
  wget https://github.com/CGAL/cgal/releases/download/v6.1-beta1/CGAL-6.1-beta1.tar.xz
  tar xf "$cgal"
  cd ${cgal%.tar.xz}
  cmake -DCMAKE_INSTALL_PREFIX=../build -DCGAL_DIR=${cgal%.tar.xz}/ -DCMAKE_BUILD_TYPE=Release -DWITH_CGAL_ImageIO=OFF -DWITH_CGAL_Qt6=OFF . 
  make install
  cd ..
fi


g++ -std=c++17 -O3 cgal.cpp -o cgal -I./build/include -L./build/lib -lmpfr -lgmp

# Create a WKT polygon from JS
node save_wkt.cjs $(cat trace.txt) poly.wkt
./cgal poly.wkt
