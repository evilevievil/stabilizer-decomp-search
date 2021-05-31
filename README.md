# stabilizer-decomp-search
Stabilizer Decomposition Searcher


##################################################
## Testing:
tons_CH_gate_test.m


##################################################
## Commands:
# allocate new state: 
CH_state(#qubits)

# initialize allocated state: 
stab.CH_init('type of initialization')

# apply gate
stab.CH_gate('type of gate', target_qubits)

# compute basis stab state inner product
CH_basis_inner_product(basis state as integer,stab)
