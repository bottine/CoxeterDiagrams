Dependency of [VinbergsAlgorithmNF](github.com/bottine/VinbergsAlgorithmNF) dealing with (hyperbolic) Coxeter diagrams.
The reference is Guglielmetti's thesis and [CoxIter](https://github.com/rgugliel/CoxIter).

Supports:

* Hyperbolic diagrams of up to 256 vertices
* of rank up to 20
* Checking finite volume, cocompactness, and computing the f-vector

Todo:

* Optionally accept Gram matrix as input.
* Check degeneracy of the diagram
* Remove the (arbitrary) limit on number of vertices and rank.
* Check whether diagram is degenerate of given rank
* Allow computing growth series and arithmeticity, etc (i.e. feature parity with CoxIter)
