* in `extend!_all_extensions`, make the code that constructs the degree sequence and associated data into a function, for clarity.
* Make the vectors `deg_1_vertices` and others into a custom "small set" type and use it everywhere  to reduce allocations etc
* Drop more bad diagrams in `extend!_all_extensions` if possible
* Function to compute all spherical subdiagrams and all affine subdiagrams
* Function to compute the f-vector

* For the following procedures:
    * compute all affines of rank from n to m
    * compute the f-vector
    * compute all direct spherical extensions
    * compute all direct affine extensions
    * check that all affine extend to at least one of rank d-1
  Follow the same pattern as the enumeration of spherical: that is, make it into a proper iterator.
  Then, this allows short-circuiting finite_volume and cocompactness tests  
  It also allows computing the f-vector by simply iterating over all spherical (no matter the rank) and incrementing the f-vector's coordinates as necessary (+ computing the affines)
  
* More unit testing:
    * Test all functions on the one specific example already present
* Write `is_degenerate()` to check if the diagram is degenerate in the given rank

# Done

* Graph Isomorphism (the best approach is probably to use [Nauty.jl](https://github.com/bovine3dom/Nauty.jl/))
  (Done using LightGraphs's experimental graph isomorphism check)

