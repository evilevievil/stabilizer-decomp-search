%%
%%
%% globals
%% CAVEAT!!!! Clear workspace before each test run!
len = 2;
vec_len = 2.^len;
num_gates = 100;
bit_H = 2.^(-0.5)*[1,1;1,-1];
bit_S = [1,0;0,1i];
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
s = CH_state(len);
s.CH_init('zero');
%% init zero state vector
state_vector = zeros(vec_len,1);
% start with uniform superposition state
for i = 1:len
    s.CH_gate('HL',i);
end

for i = 1:vec_len
    state_vector(i) = sqrt(1/vec_len); 
end

gate_record = zeros(num_gates,3);

s_conj = CH_state(len);
    s_conj.transpose(s);
    s_conj.pp_CH('ch');
    s_conj.pp_CH('basis');

%% apply gates to construct s
for i = 1:num_gates
    gate_choice = randi(2,1,1);
    if gate_choice == 1  %SLs
        bit_choice = randi(len,1,1);
        gate = S_array(:,:,bit_choice);
        s.CH_gate('SL',bit_choice);
        gate_record(i,:) = [gate_choice,bit_choice(1),-1];
    elseif gate_choice == 2 %CXL
        bit_choice = randperm(len,2); % use 1st as control and 2nd as result
        gate = CX_array(:,:,bit_choice(1),bit_choice(2));
        s.CH_gate('CXL',bit_choice);
        gate_record(i,:) = [gate_choice,bit_choice(1),bit_choice(2)];
    elseif gate_choice == 3 %CZL
        bit_choice = randperm(len,2); % use 1st as control and 2nd as result
        gate = CZ_array(:,:,bit_choice(1),bit_choice(2));
        s.CH_gate('CZL',bit_choice);
        gate_record(i,:) = [gate_choice,bit_choice(1),bit_choice(2)];
    else
        fprintf('error invalid gate choice.\n');
    end
    fprintf('%dth test\n',i);
    s_conj = CH_state(len);
    s_conj.transpose(s);
    s_conj.pp_CH('ch');
    s_conj.pp_CH('basis');
end

%%
s.pp_CH('ch');
s.pp_CH('basis');
s_conj = CH_state(len);
s_conj.transpose(s);
s_conj.pp_CH('ch');
s_conj.pp_CH('basis');

%% apply gates to transpose
for i = 1:num_gates
    gate_choice = gate_record(i,1);
    bit_choice_1 = gate_record(i,2);
    bit_choice_2 = gate_record(i,3);
    if gate_choice == 1  %SLs
        s_conj.CH_gate('SL',bit_choice_1);
    elseif gate_choice == 2 %CXL
        s_conj.CH_gate('CXL',[bit_choice_1,bit_choice_2]);
    elseif gate_choice == 3 %CZL
        s_conj.CH_gate('CZL',[bit_choice_1,bit_choice_2]);
    else
        fprintf('error invalid gate choice.\n');
    end
    fprintf('%dth test\n',i);
    %s_conj.pp_CH('ch');
    s_conj_state_vec = CH2basis(s_conj);
    %disp(state_vector); disp(s_conj_state_vec); 
end

assert(approx_equal(state_vector,s_conj_state_vec,0.000000001));
fprintf('conjugate test passed!\n');

