## The example shows:
- How to use the *raw* interface of the bindings.
- How to use `fetch_many()` and `push_many()`.
- Simple profiling for `fetch()`and `push()`.

## Notes:
- To run the example `bash run_case.sh`.

> Nt (time steps), Ni (i coords) and Nj (j coords) can be changed in the script.

- Raw interface could be around  **3x-5x faster** than Pythonic one for *fetch/push*. This depends much on the python extra code overhead *VS* samplers/communication time.
- Raw interface does not perform any type checking.
- `push_many()`/`fetch_many()` can be 1 or 2 orders of magnitude faster than python looping using `fetch()`/`push()`at the expense of having to store points and values first in a numpy array.

