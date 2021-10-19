# Efficient Network Reliability Estimation for Client-Server Model

This repository includes the codes to reproduce the results of experiments in the paper "Efficient Network Reliability Evaluation for Client-Server Model" published in the proceedings of [2021 IEEE Global Communications Conference (GLOBECOM 2021)](https://globecom2021.ieee-globecom.org/).

## Requirements

All codes are written in C++11 language. Before building our code, you must place the header files of [TdZdd](https://github.com/kunisura/TdZdd) into `src/` directory. More specifically, all header files of [TdZdd/include/tdzdd](https://github.com/kunisura/TdZdd/tree/master/include/tdzdd) must be placed on `src/tdzdd/` directory; e.g. `src/tdzdd/DdEval.hpp`, `src/tdzdd/DdSpec.hpp`, etc.

After installing TdZdd, if your environment has CMake version >=3.8, you can build all codes with the following command:

```shell
(moving to src/ directory)
mkdir release
cd release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

After running this, all binaries are generated at the `release/` directory. Alternatively, since our code uses no external library (other than TdZdd), you can build each binary separately by your compiler command like:

```shell
g++ -std=c++11 -O3 -DNDEBUG -o main main.cpp mylib/graph.cpp
g++ -std=c++11 -O3 -DNDEBUG -o main_single main_single.cpp mylib/graph.cpp
g++ -std=c++11 -O3 -DNDEBUG -o tdzdd tdzdd.cpp mylib/graph.cpp
g++ -std=c++11 -O3 -DNDEBUG -o tdzdd_single tdzdd_single.cpp mylib/graph.cpp
```

### Verified environments

We verified that the building process of our codes and the commands presented below worked fine in the following macOS and Linux environments:

- macOS Big Sur 11.2.1 + Apple clang 12.0.0
- CentOS 7.9 + gcc 4.8.5
- CentOS 7.3 + gcc 9.1.0

## How to reproduce experimental results

After building our code, the following 4 binaries are generated: `main`, `main_single`, `tdzdd`, and `tdzdd_single`. `main` and `main_single` implements our proposed method, while `tdzdd` and `tdzdd_single` implements the existing (HH) method. The suffix `_single` means that it solves the classic $k$-NR problem instead of the new $k$-NR+ problem described in our paper.

All data used in our experiments are in `data.tar.gz`. After extracting, it includes the data of the following graphs:
- Grid graphs: `grid7x14.txt` (Grid-7x14), `grid7x28.txt` (Grid-7x28), `grid7x42.txt` (Grid-7x42), `grid8x8.txt` (Grid-8x8), `grid10x10.txt` (Grid-10x10), `grid12x12.txt` (Grid-12x12)
- [Topology Zoo](http://www.topology-zoo.org/index.html): `0146-real-Interoute.edgelist.txt` (Interroute), `0186-real-TataNld.edgelist.txt` (TATA), `0895-real-Kdl.edgelist.txt` (Kentucky Datalink)
- [Rocketfuel](https://research.cs.washington.edu/networking/rocketfuel/): `1221.dat.txt` (Rocketfuel-1221), `1755.dat.txt` (Rocketfuel-1755), `3257.dat.txt` (Rocketfuel-3257), `3967.dat.txt` (Rocketfuel-3967), `6461.dat.txt` (Rocketfuel-6461)

For each graph, the following data are included:

- `<graphname>`: the edge list of the graph.
- `<graphname>.prob`: the working probability ($p_i$) of each edge ($e_i$) in the same order as the edge list. As described in our paper, each $p_i$ is chosen uniformly at random from $[0.9,0.95]$.
- `<graphname>.bc<k>.src`: the list of `<k>` vertices with higher betweenness centrality.

_Note:_ Since our primary purpose is to confirm the computational cost of our and existing methods, we preliminarily remove degree 1 vertices from these graph data as described in Section II of our paper.

### Experiments for $k$-NR+

To check the computational cost of our method for solving $k$-NR+, run:

```shell
./main [graph_file] [probability_file] [source_vertices_file] [order_file]
```

This command executes our proposed method with graphs, probabilities and source vertices specified by `[graph_file]`, `[probability_file]` and `[source_vertices_file]`, respectively. `[order_file]` specifies the edge ordering used in our method. When using the data of this repository, `[order_file]` will be the same as `[graph_file]`.

The exsiting method can also be executed as follows:

```shell
./tdzdd [graph_file] [probability_file] [source_vertices_file] [order_file]
```

_Running example:_ To execute our proposed method for the Grid-7x14 graph with source the $k=$5 vertices with higher betweeness centrality, run:

```shell
./main data/grid7x14.txt data/grid7x14.txt.prob data/grid7x14.txt.bc5.src data/grid7x14.txt
```

### Experiments for $k$-NR

To check the computational cost of our method for solving $k$-NR, run:

```shell
./main_single [graph_file] [probability_file] [terminal_vertices_file] [order_file]
```

This command executes our proposed method with graphs, probabilities and terminal vertices specified by `[graph_file]`, `[probability_file]` and `[terminal_vertices_file]`, respectively. `[order_file]` specifies the edge ordering used in our method. When using the data of this repository, `[order_file]` will be the same as `[graph_file]`. When running our method for $k$-NR, the first $k-1$ vertices are regarded as source vertices and the last one vertex is regarded as an destination.

The exsiting method can also be executed as follows:

```shell
./tdzdd_single [graph_file] [probability_file] [terminal_vertices_file] [order_file]
```

_Running example:_ To execute our proposed method for the Grid-7x14 graph with terminal the $k=$6 vertices with higher betweeness centrality, run:

```shell
./main_single data/grid7x14.txt data/grid7x14.txt.prob data/grid7x14.txt.bc6.src data/grid7x14.txt
```

## License

This software is released under the NTT license, see `LICENSE.txt`.
