#! /bin/sh

# sudo apt-get install libcgal-dev libcgal-qt5-dev libgmp-dev libmpfr-dev

# Compile (typical):
# g++ -std=c++17 -O2 cgal.cpp -o cgal $(pkg-config --cflags --libs cgal) -lgmp -lmpfr

# Or without pkg-config (may vary by distro):
g++ -std=c++17 -O3 cgal.cpp -o cgal -lmpfr -lgmp

# Create a WKT polygon from JS
node save_wkt.cjs $(cat trace.txt) poly.wkt
./cgal poly.wkt
