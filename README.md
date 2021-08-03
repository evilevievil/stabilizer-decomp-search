## stabilizer-decomp-search
Stabilizer Decomposition Searcher

# Testing:
Unit tests are stored under 'unit_tests' folder. Note that 'tons_stab_decomp_proj_test' is currently uncompatible with the memoization optimization. Functions are tested using state vector simulations.
1. tons_CH_gate_test.m (1000 random gates test)
2. tons_CH_conjugate_test.m (random gates test for conjugate; equal superposition initial state; tests without HL gate)
3. tons_CH_conjugate_HL_test.m (random gates test for conjugate; tests with HL gate)
4. tons_CH_pauli_proj_test.m (random pauli projector test)
5. tons_stab_decomp_proj_test.m (random stabilizer decomposition projection + pauli projector update test)
6. CH_test.m (sanity/mock tests; please configure 'const.m' to max 8 qubits when testing)

# Configuration:
1. max number of qubits: in 'const.m', update init_max_qubits = 8 or 16 or 32

# Commands:
1. allocate new state: 
```c
stab = CH_state(#qubits) % #qubits <= init_max_qubits
```
2. initialize allocated state: 
```c
stab.CH_init('type_of_initialization')
% type_of_initialization = {'zero' (|00..0>),'rand' (random stab state)}
```
3. apply gate: 
```c
stab.CH_gate('type_of_gate', target_qubits)
% type_of_gate = {'CXL','CXR','CZL','CZR','HL','SL','SR'}
% target_qubits = {[control_bit,target_bit],[trarget_bit,dont_care]}
```
4. apply pauli projector: \
*** CAVEAT: the function normalizes the result for you; no need to renormalize 
```c
stab.CH_pauli_proj(is_neg,x_bit,z_bit)
% is_neg = positive->0 ; negative->1
% x_bit[i] = 0->I ; 1->X
% z_bit[i] = 0->I ; 1->Z
```
5. compute basis stab state inner product: \
*** CAVEAT: Please follow reverse bit string convention for basis state \
*** eg. basis state'1011' should be inputted as '1101' in base 10 
```c
CH_basis_inner_product(basis state as integer,stab)
```
6. compute stab state stab state inner product: \
*** CAVEAT: the function conjugates |stab_state1> for you: please input |stab_state1>, |stab_state2> 
```c
CH_CH_inner_product(stab_state1,stab_state2)
```
7. pretty print state: 
```c
%% print CH-form 
state.pp_CH('ch') 
%% print conjugate state's CH-form 
state.pp_CH('conj') 
%% print as state vector, uses regular (not reversed) bit string convention for your convenience 
state.pp_CH('basis') 
```
8. compute projection onto a stab decomp: \
*** CAVEAT: base case for its memoized version 'CH_decomp_project_memoize' 
```c
[projection,memoize_G,memoize_target_inner_prod] = CH_decomp_project(target_state_vec,stab_decomp,#qubit,decomp_len)
```
9. search for stabilizer decomposition given target_state and target_rank:
```c
fixed_rank_stab_decomp_search(target_state,#qubits,target_rank,init_temp_inverse,final_temp_inverse,max_SA_steps,rand_walk_steps)
```
10. make state vector for magic state:
```c
state_vec = magic_state_vec(type_of_magic, #qubits)
% type_of_magic = {'T', 'H', 'catT'}
```
