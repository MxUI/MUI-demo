## The example shows
- How to use `create_unifaces()` to create a group of same type interfaces.
- How to use `set_data_types_unifaces()` to assign types to data tags all at once.
-  Define multiple interfaces of different types between different codes.

## Topology interface map of the case
```bash
                        (N=1 ranks)
                         [code1-1]
                           /
                       {I1-2D}
 (N=1 ranks)             /
  [code2-3]--{I3-3D}--[code0-0] (N=2 ranks)
                         \
                       {I2-2D}
                           \     
                        [code1-2]--{I4-1D}--[code3-4] 
                       (N=2 ranks)         (N=1 ranks)
    [*] Legend:
         - {Ii-kD}: k-dimensional interface with label i
         - [codei-j]: Instances of a code of type i associated with domain j (N processes spawned)
```
## Notes
- A single interface can only be shared by a pair of codes.
- All codes has to call `create_unifaces()` for each type of interface used among all processes, even if a code itself do not use that particular interface. This is due to the need of an internal collective MPI call on the **MPI_COMM_WORLD** communicator.
