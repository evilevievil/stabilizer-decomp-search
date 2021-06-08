%%
%%
%% globals
len = 2;
vec_len = 2.^len;
num_gates = 1000;
bit_H = 2.^(-0.5)*[1,1;1,-1];
bit_S = [1,0;0,1i];
bit_ST = [1,0;0,1i];
bit_X = [0,1;1,0];
bit_Z = [1,0;0,-1];
bit_I = [1,0;0,1];
bit_I_0 = [1,0;0,0];
bit_I_1 = [0,0;0,1];
rng('default');
%% make gates by tensor product
S_array = zeros(vec_len,vec_len,len);
H_array = zeros(vec_len,vec_len,len);
CX_array = zeros(vec_len,vec_len,len,len);
CZ_array = zeros(vec_len,vec_len,len,len);
for i = 1:len
    S_array(:,:,i) = kron(kron(tensor_exp(bit_I,i-1),bit_S),tensor_exp(bit_I,len-i));
    H_array(:,:,i) = kron(kron(tensor_exp(bit_I,i-1),bit_H),tensor_exp(bit_I,len-i));
end

for i = 1:len
    for j = 1:len
        if i==j 
            continue;
        elseif i<j 
            CX_ij_I = kron(kron(kron(kron(tensor_exp(bit_I,i-1),bit_I_0),tensor_exp(bit_I,j-i-1)),bit_I),tensor_exp(bit_I,len-j));
            CX_ij_X = kron(kron(kron(kron(tensor_exp(bit_I,i-1),bit_I_1),tensor_exp(bit_I,j-i-1)),bit_X),tensor_exp(bit_I,len-j));
            CZ_ij_I = kron(kron(kron(kron(tensor_exp(bit_I,i-1),bit_I_0),tensor_exp(bit_I,j-i-1)),bit_I),tensor_exp(bit_I,len-j));
            CZ_ij_Z = kron(kron(kron(kron(tensor_exp(bit_I,i-1),bit_I_1),tensor_exp(bit_I,j-i-1)),bit_Z),tensor_exp(bit_I,len-j));
            CX_array(:,:,i,j) = CX_ij_I + CX_ij_X;
            CZ_array(:,:,i,j) = CZ_ij_I + CZ_ij_Z;
        else
            CX_ij_I = kron(kron(kron(kron(tensor_exp(bit_I,j-1),bit_I),tensor_exp(bit_I,i-j-1)),bit_I_0),tensor_exp(bit_I,len-i));
            CX_ij_X = kron(kron(kron(kron(tensor_exp(bit_I,j-1),bit_X),tensor_exp(bit_I,i-j-1)),bit_I_1),tensor_exp(bit_I,len-i));
            CZ_ij_I = kron(kron(kron(kron(tensor_exp(bit_I,j-1),bit_I),tensor_exp(bit_I,i-j-1)),bit_I_0),tensor_exp(bit_I,len-i));
            CZ_ij_Z = kron(kron(kron(kron(tensor_exp(bit_I,j-1),bit_Z),tensor_exp(bit_I,i-j-1)),bit_I_1),tensor_exp(bit_I,len-i));
            CX_array(:,:,i,j) = CX_ij_I + CX_ij_X;
            CZ_array(:,:,i,j) = CZ_ij_I + CZ_ij_Z;
        end
    end
end



%% init ch zero state
s1 = CH_state(len);
s1.CH_init('zero');
s2 = CH_state(len);
s2.CH_init('zero');
%% init zero state vector
state_vector1 = zeros(vec_len,1);
state_vector1(1) = 1; % zero state
state_vector2 = zeros(vec_len,1);
state_vector2(1) = 1; % zero state
state_vector3 = zeros(vec_len,1);
state_vector3(1) = 1; % zero state


gate_record = zeros(num_gates,3);

%% test left gates
for i = 1:num_gates
    fprintf('%dth ch ch inner product test!\n',i);
    gate_choice1 = randi(6,1,1);
    gate_choice2 = randi(6,1,1);
    if gate_choice1 == 1  %SLs
        bit_choice = randi(len,1,1);
        gate1 = S_array(:,:,bit_choice);
        s1.CH_gate('SL',bit_choice);
        fprintf('gate1:SL on bit: %d\n',bit_choice(1));
    elseif gate_choice1 == 2  %HL
        bit_choice = randi(len,1,1);
        gate1 = H_array(:,:,bit_choice);
        s1.CH_gate('HL',bit_choice);
        fprintf('gate1:HL on bit: %d\n',bit_choice(1));
    elseif gate_choice1 == 3 || gate_choice1 == 4 %CXL
        bit_choice = randperm(len,2); % use 1st as control and 2nd as result
        gate1 = CX_array(:,:,bit_choice(1),bit_choice(2));
        s1.CH_gate('CXL',bit_choice);
        fprintf('gate1:CXL on bits: %d, %d\n',bit_choice(1),bit_choice(2));
    elseif gate_choice1 == 5 || gate_choice1 == 6 %CZL
        bit_choice = randperm(len,2); % use 1st as control and 2nd as result
        gate1 = CZ_array(:,:,bit_choice(1),bit_choice(2));
        s1.CH_gate('CZL',bit_choice);
        fprintf('gate1:CZL on bits: %d, %d\n',bit_choice(1),bit_choice(2));
    else
        fprintf('error invalid gate choice.\n');
    end

    if gate_choice2 == 1  %SLs
        bit_choice = randi(len,1,1);
        gate2 = S_array(:,:,bit_choice);
        s2.CH_gate('SL',bit_choice);
        fprintf('gate2:SL on bit: %d\n',bit_choice(1));
    elseif gate_choice2 == 2  %HL
        bit_choice = randi(len,1,1);
        gate2 = H_array(:,:,bit_choice);
        s2.CH_gate('HL',bit_choice);
        fprintf('gate2:HL on bit: %d\n',bit_choice(1));
    elseif gate_choice2 == 3 || gate_choice2 == 4 %CXL
        bit_choice = randperm(len,2); % use 1st as control and 2nd as result
        gate2 = CX_array(:,:,bit_choice(1),bit_choice(2));
        s2.CH_gate('CXL',bit_choice);
        fprintf('gate2:CXL on bits: %d, %d\n',bit_choice(1),bit_choice(2));
    elseif gate_choice2 == 5 || gate_choice2 == 6 %CZL
        bit_choice = randperm(len,2); % use 1st as control and 2nd as result
        gate2 = CZ_array(:,:,bit_choice(1),bit_choice(2));
        s2.CH_gate('CZL',bit_choice);
        fprintf('gate2:CZL on bits: %d, %d\n',bit_choice(1),bit_choice(2));
    else
        fprintf('error invalid gate choice.\n');
    end

    state_vector1 = gate1 * state_vector1;
    state_vector2 = gate2 * state_vector2;
    state_vector_prod = dot(state_vector1,state_vector2);
    ch_ch_prod = CH_CH_inner_product(s1,s2);

    %disp(state_vector); disp(s_state_vec);
    %s.pp_CH('ch');
    assert(approx_equal(state_vector1,CH2basis(s1),0.000000001));
    assert(approx_equal(state_vector2,CH2basis(s2),0.000000001));
    disp(state_vector_prod); disp(ch_ch_prod);
    disp(state_vector1); disp(state_vector2);
    assert(approx_equal(CH_CH_inner_product(s1,s1),1,0.000000001));
    assert(approx_equal(CH_CH_inner_product(s2,s2),1,0.000000001));
    assert(approx_equal(state_vector_prod,ch_ch_prod,0.000000001)); %% may need +- to account for rounding error...
    fprintf('%dth ch ch inner product test passed!\n',i);
end
