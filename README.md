# pph-cpp

Developer: Samir Chowdhury

Date: March 26, 2019

Paper: Chowdhury, S. and Mémoli, F., Persistent Path Homology of Directed Networks. SODA 2018.

This is a C++11 implementation of the Persistent Path Homology (PPH) package described in the aforementioned paper. 

Tested using the following compiler:
Apple LLVM version 10.0.0 (clang-1000.10.44.4)

### Compilation
Compile by running the following in terminal:

clang++ -std=c++11 computePPH.cpp -o computePPH

For optimization, use the following:

clang++ -std=c++11 -Ofast -march=native computePPH.cpp -o computePPH


### Usage

The .txt files provided in the tests folder are examples of input files to pph-cpp. 
The first line in the input file should be the number of nodes. In the following
lines, each row has three entries: the first two describe the edge, and the 
third entry is the corresponding edge weight.

After compilation, use the code as follows (in a terminal window):

./computePPH

*Type in a file name (e.g. 'cycleNet10.txt')*

cycleNet10.txt

*done! pers1 is *

*1 5 *

*time elapsed 21 milliseconds*


The output is saved in csv format as res_cycleNet10.txt.
