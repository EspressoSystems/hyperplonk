KZG based multilinear polynomial commitment
-----

# Compiling features:
- `parallel`: use multi-threading when possible.
- `print-trace`: print out user friendly information about the running time for each micro component.
- `group-switched`: witch the groups of G1 and G2 as in a normal KZG. i.e., commitments are over G2 and openings are over G1 when this feature is turned on.