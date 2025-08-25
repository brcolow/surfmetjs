#! /bin/sh

# sudo apt-get install libcgal-dev libcgal-qt5-dev libgmp-dev libmpfr-dev hotspot

cgal="CGAL-6.1-beta1.tar.xz"
if ! [ -f "$cgal" ]; then
  wget https://github.com/CGAL/cgal/releases/download/v6.1-beta1/CGAL-6.1-beta1.tar.xz
  tar xf "$cgal"
  cd ${cgal%.tar.xz}
  cmake -DCMAKE_INSTALL_PREFIX=../build -DCGAL_DIR=${cgal%.tar.xz}/ -DCMAKE_BUILD_TYPE=Release -DWITH_CGAL_ImageIO=OFF -DWITH_CGAL_Qt6=OFF . 
  make install
  cd ..
fi


g++ -std=c++17 -O3 -g -fno-omit-frame-pointer cgal.cpp -I./build/include -L./build/lib -lmpfr -lgmp -o cgal

# Create a WKT polygon from nist trace
node save_wkt.cjs ../nist/cir2d10.ds poly.wkt
# sudo ./perf-6.9.0/tools/perf/perf record -o ./cgal.perf.data -F 99 -g -- ./cgal poly.wkt
# hotspot perf.data
./cgal poly.wkt
