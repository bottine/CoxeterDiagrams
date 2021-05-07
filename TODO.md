* in `extend!_all_extensions`, make the code that constructs the degree sequence and associated data into a function, for clarity.
* Make the vectors `deg_1_vertices` and others into a custom "small set" type and use it everywhere  to reduce allocations etc
* Drop more bad diagrams in `extend!_all_extensions` if possible
* make `extend!` into a non-recursive function
* make `all_affine_extend_well` into a non-recursive iterating function and compare to the current naïve implementation

  
* More unit testing:
    * Test all functions on the one specific example already present
* Write `is_degenerate()` to check if the diagram is degenerate in the given rank

# Done

* Graph Isomorphism (the best approach is probably to use [Nauty.jl](https://github.com/bovine3dom/Nauty.jl/))
  (Done using LightGraphs's experimental graph isomorphism check)

