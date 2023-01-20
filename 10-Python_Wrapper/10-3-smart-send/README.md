[MUI minimum required: (VERSION 1.2+)]
## The example shows:
- How to use *smart send* capabilities of MUI with different span shapes and timeouts.
- How to use `mpi4py` in combination with `mui4py` module.

## Topology
```
    [*] Decomposition of the domain [-1,-1] x [1,1] for receiver processes
        x: Fetched points at the center of quadrants
        Ri: Ranks

     (-1,1) |=====|=====| (1,1)
            | R2  | R3  |
            |  x  |  x  |
            |     |     |
            |=====|=====|
            | R0  | R1  |
            |  x  |  x  |
            |     |     |
    (-1,-1) |=====|=====| (1,-1)

    [*] Decomposition of the domain [-1,-1] x [1,1] for sender processes
        x: pushed points at the center of quadrants
        Ri: Ranks


     (-1,1) |===========| (1,1)
            |           |
            |  x     x  |
            |           |
            |     R0    |
            |           |
            |  x     x  |
            |           |
    (-1,-1) |===========| (1,-1)


     Eeach loop iteration (step from 1 to 3) sender uses a different span:
         s = 1 -> Box intersecting R0 and R2
         s = 2 -> Sphere intersecting R1 and R3
         s = 3 -> Point intersecting R3
         s = 4 -> Smart sending is disabled due to timeout


```

## Notes:
- To run the example `bash run.sh`.
- Empty fetched values indicate the value is not available for the queried point in the rank.
