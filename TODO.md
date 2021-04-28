* in `extend!_all_extensions`, make the code that constructs the degree sequence and associated data into a function, for clarity.
* Make the vectors `deg_1_vertices` and others into a custom "small set" type and use it everywhere  to reduce allocations etc
* Drop more bad diagrams in `extend!_all_extensions` if possible
* Function to compute all spherical subdiagrams and all affine subdiagrams
* Function to compute the f-vector

* Make `all_spherical_of_rank` into an iterator using `@resumable` (same for `all_affine_of_rank` and `all_spherical_direct_extensions` and `all_affine_direct_extensions`) so that:
    * Iterating allocates less
    * Can short-circuit the iteration as soon as needed
* More unit testing:
    * Test all functions on the one specific example already present
* Write `is_degenerate()` to check if the diagram is degenerate in the given rank

# Done

* Graph Isomorphism (the best approach is probably to use [Nauty.jl](https://github.com/bovine3dom/Nauty.jl/))
  (Done using LightGraphs's experimental graph isomorphism check)

