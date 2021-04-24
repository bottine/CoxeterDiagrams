* Make `all_spherical_of_rank` into an iterator using `@resumable` (same for `all_affine_of_rank` and `all_spherical_direct_extensions` and `all_affine_direct_extensions`) so that:
    * Iterating allocates less
    * Can short-circuit the iteration as soon as needed
* More unit testing:
    * Test all functions on the one specific example already present

# Done

* Graph Isomorphism (the best approach is probably to use [Nauty.jl](https://github.com/bovine3dom/Nauty.jl/))
  (Done using LightGraphs's experimental graph isomorphism check)

