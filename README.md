# LDSieve
##  Implementation of lattice sieving with List Decoding
**version 0.9** -- *Release date: September 2015*

**Author:** Leo Ducas

LDSieve is open-source software distributed under the terms of the GNU
General Public License. See the file LICENSE for complete details on the licensing. LDSieve is a prototype implementation of the Lattice Sieve Algorithm
described in the paper

["New directions in nearest neighbor searching with applications to lattice sieving",  Anja Becker and Léo Ducas and Nicolas Gama and Thijs Laarhoven]

The code is released as an effort to stimulate new research and efficient  implementations in the concerned fields, namely Lattice reduction algorithm,
and Nearest Neighbor Search algorithms. The codes comes with no warranty. Additionally, we provide the benchmarking procedure and the instances used in the paper for reproducibility.

### Requirements 

The code require NTL to be installed in your home directory (if installed somewhere else, update the Makefile).

The Benchmark script requires python.


### Installation

Make.

### Running Benchmark

In the Benchmark/ subdirectory, run python Benchmark.py .44 .44 > 44.bench

