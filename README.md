## stabilizer-decomp-search
Stabilizer Decomposition Searcher

# Testing:
1. tons_CH_gate_test.m (1000 random gates test)
2. CH_test.m (sanity/mock tests)

# Commands:
1. allocate new state: \
CH_state(#qubits)

2. initialize allocated state: \
stab.CH_init('type of initialization')

3. apply gate: \
stab.CH_gate('type of gate', target_qubits)

4. compute basis stab state inner product: \
CH_basis_inner_product(basis state as integer,stab)
