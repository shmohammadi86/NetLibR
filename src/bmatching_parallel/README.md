This repo contains a school project of mine: the task was to write, fine-tune and evaluate an efficient parallel implementation of an *b*-suitors algorithm [(Khan et al., 2016)](https://www.cs.purdue.edu/homes/apothen/Papers/bMatching-SISC-2016.pdf) for solving the problem of finding a maximal *b*-matching of a weighted undirected graph. The problem of *b*-matching generalizes the graph matching problem: the goal is to match each vertex with at most *b* other vertices in a way that maximizes their weights. This problem has some interesting applications in privacy-aware machine learning (cf. [adaptive anonymity](https://papers.nips.cc/paper/4858-adaptive-anonymity-via-b-matching.pdf)).

## Design outline

A C++14 implementation of *b*-suitors algorithm was prepared that can be provided with an arbitrary pre-defined (maximal) number of threads and an arbitrary assignment of *b* to each vertex. In its main loop, the algorith iterates over vertices trying to find the best matches, until a vertex finds *b* matches or there is no eligible candidate left. This loop was implemented to run in parallel for each vertex. Threads are joined at the end of each full iteration over the graph. Additionally, exclusive access to the matching queue was ensured for each thread when finalizing the matching.

The scheduling regime employed was a thread pool with fixed number of threads and FIFO task queue.

The implementation was evaluated on datasets from [Stanford Large Network Dataset Collection](http://snap.stanford.edu/data/).

## Quickstart

```bash
cmake .
./adorate thread-count inputfile b-limit
```

`inputfile` is a path to file containing the graph (each line specifies two vertices ids and a weight). ``b-limit`` is the maximal *b*. `blimit.cpp` implements an assignment of *b* to each vertex.
