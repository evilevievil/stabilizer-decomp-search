%%
%%
%% globals
len = 2;
vec_len = 2.^len;
decomp_len = 2;
num_iter = 100;
a = magic_state_vec('T',len);

% reformat amplitude array to reverse bit string convention
reverse_formatted_a = zeros(vec_len,1);
for i = 1:vec_len
    %% matlab 1 indexing!! :(
    %% reverse bit string due right most convention
    diff_len = len - strlength(dec2bin(i-1));
    i_strrev = reverse(dec2bin(i-1));
    for j = 1: diff_len
        i_strrev = append(i_strrev,'0');
    end
    reverse_formatted_a(i,1) = a(bin2dec(i_strrev)+1,1);
end

%% init stab_decomp
for i= 1:decomp_len
    stab_decomp(i) = CH_state(len);
    stab_decomp(i).CH_init('rand');
end

%% stab decomp proj test
for i = 1:num_iter
    fprintf('%d stab decomp proj test\n', i);
    pauli_update(stab_decomp,decomp_len);
    ch_norm = CH_decomp_project(reverse_formatted_a,stab_decomp,len,decomp_len);
    state_vec_norm = state_vector_proj(a,stab_decomp,len,decomp_len);
    disp(state_vec_norm);
    disp(ch_norm);
    assert(approx_equal(state_vec_norm,ch_norm,0.000000001));
    fprintf('%d stab decomp proj test passed!\n', i);
end

%% 
function proj_norm = state_vector_proj(a,stab_decomp, len, decomp_len)
    vec_len = 2.^len;
    state_vec_decomp = zeros(vec_len,decomp_len);
    for i = 1:decomp_len
        state_vec_decomp(1:vec_len,i) = CH2basis(stab_decomp(i));
    end
    orth_decomp = orth(state_vec_decomp);
    state_vec_proj = orth_decomp * orth_decomp';
    %disp(state_vec_proj);
    %disp(state_vec_proj*state_vec_proj);
    assert_result = approx_equal(state_vec_proj,state_vec_proj*state_vec_proj,0.000000001);
    for i = 1:decomp_len
        assert(assert_result(i));
    end
    proj_norm = a' * state_vec_proj * a;
end

%% 
function pauli_update(stab_decomp,decomp_len)
    state_choice = randi(decomp_len,1,1);
    stab_decomp(state_choice) = stab_decomp(state_choice).CH_pauli_proj('rand',-1);
end
