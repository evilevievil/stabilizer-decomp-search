d%%
%%
%% globals
len = 8;
vec_len = 2.^len;
decomp_len = 6;
num_iter = 100;
a = magic_state_vec('T',len);
%a = [0.5;0.5;0.5;-0.5i];
%rng('default');
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
    new_stab_decomp = pauli_update(stab_decomp,decomp_len);
    ch_norm = CH_decomp_project(reverse_formatted_a,new_stab_decomp,len,decomp_len);
    if ch_norm < 0
       i = i-1;
       continue;
    end
    stab_decomp = new_stab_decomp;
    state_vec_norm = state_vector_proj(a,stab_decomp,len,decomp_len);
    disp(state_vec_norm);
    disp(ch_norm);
    for j = 1:decomp_len
        state_vec_decomp(1:vec_len,j) = CH2basis(stab_decomp(j));
    end
    %disp(state_vec_decomp)
    assert(approx_equal(state_vec_norm,ch_norm,0.000000001));
    fprintf('%d stab decomp proj test passed!\n', i);
end

%% Test with Gram matrix state vec projection
function proj_norm = state_vector_proj1(a,stab_decomp, len, decomp_len)
    vec_len = 2.^len;
    a_stab_array = zeros(1,decomp_len);
    state_vec_decomp = zeros(vec_len,decomp_len);
    for i = 1:decomp_len
        state_vec_decomp(1:vec_len,i) = CH2basis(stab_decomp(i));
    end
    proj_norm = double(0);
    
    % compute G array
    G = zeros(decomp_len);
    for i = 1:decomp_len
        % todo: just use conjugate for lower triangular
        for j=1:decomp_len
            G(i,j) = dot(state_vec_decomp(:,i),state_vec_decomp(:,j));
        end
    end
    
     if det(G) ~= 0
        G_inv = inv(G);
        % compute basis array
        for i = 1:decomp_len
            for j = 1:vec_len
                % careful! Matlab 1 indexing :(
                a_stab_array(i) = a_stab_array(i)  + conj(a(j)) * state_vec_decomp(j,i);
                %disp(a_stab_array);
            end
        end
        
        %disp(a_stab_array);
        %disp(G);

        % compute projection
        for i = 1:decomp_len
            for j = 1:decomp_len
                proj_norm = proj_norm + a_stab_array(:,i) * conj(a_stab_array(:,j)) * G_inv(i,j);
            end
        end
    else
        % Gram matrix is singular -> find a better stabilizer decomp
        proj_norm = -1;
    end
end

%% Test with orth projection state vec projection
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
    %assert_result = approx_equal(state_vec_proj,state_vec_proj*state_vec_proj,0.000000001);
    %for i = 1:decomp_len
    %    assert(assert_result(i));
    %end
    proj_norm = a' * state_vec_proj * a;
    if rank(state_vec_proj,0.00001)<decomp_len
        proj_norm = -1;
    end
end

%% 
function new_stab_decomp = pauli_update(stab_decomp,decomp_len)
    state_choice = randi(decomp_len,1,1);
    bit_choice = randi(4,len,1);
    sign_choice = randi(2,1,1)-1;
    x_bits = const.init_uint;
    z_bits = const.init_uint;
    for i = 1:len
        if bit_choice(i,1) == 1 % I  
            continue;
        elseif bit_choice(i,1) == 2 % X
            x_bits = bitset(x_bits,i,1);
        elseif bit_choice(i,1) == 3 % Z
            z_bits = bitset(z_bits,i,1);
        else % Y
            x_bits = bitset(x_bits,i,1);
            z_bits = bitset(z_bits,i,1);
        end
    end
    stab_decomp(state_choice) = stab_decomp(state_choice).CH_pauli_proj(sign_choice,x_bits,z_bits);
    new_stab_decomp = stab_decomp;
end
