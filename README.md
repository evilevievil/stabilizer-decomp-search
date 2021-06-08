## stabilizer-decomp-search
Stabilizer Decomposition Searcher

# Testing:
1. tons_CH_gate_test.m (1000 random gates test)
2. tons_CH_conjugate_test.m (random gates test for conjugate)
3. tons_CH_pauli_proj_test.m (random pauli projector test)
4. CH_test.m (sanity/mock tests)

# Commands:
1. allocate new state: \
CH_state(#qubits)

2. initialize allocated state: \
stab.CH_init('type of initialization')

3. apply gate: \
stab.CH_gate('type of gate', target_qubits)

4. apply pauli projector \
stab.CH_pauli_proj('type of pauli projector', target_qubit)

5. compute basis stab state inner product: \
*** CAVEAT: Please follow reverse bit string convention for basis state \
*** eg. basis state'1011' should be inputted as '1101' in base 10 \
CH_basis_inner_product(basis state as integer,stab)

6. compute stab state stab state inner product: \
*** CAVEAT: the function conjugates stab_state1 for you: please input |stab_state1>, |stab_state2> \
CH_CH_inner_product(stab_state1,stab_state2)

7. pretty print state \
%% print CH-form \
state.pp_CH('ch') \
%% print as state vector, uses regular (not reversed) bit string convention for your convenience \
state.pp_CH('basis') 

