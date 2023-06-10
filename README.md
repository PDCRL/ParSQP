1. Clone the nlopt library from https://github.com/stevengj/nlopt.git
2. Replace CMakeLists.txt present in nlopt folder with our CMakeLists.txt
3. Replace slsqp.c from nlopt/src/algs/slsqp/slsqp.c with our slsqp.cpp

Steps To Build Nlopt:
1. cd nlopt
2. mkdir build
3. cd build
4. cmake ..
5. make
6. sudo make install
