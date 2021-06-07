%%
%%
%% globals
len = 3;
vec_len = 2.^len;
bit_X = [0,1;1,0];
bit_Z = [1,0;0,-1];
bit_Y = [0,-1i;1i,0];
bit_I = [1,0;0,1];
rng('default');
%% make gates by tensor product
Pauli_X_array = zeros(vec_len,vec_len,len,2);
Pauli_Z_array = zeros(vec_len,vec_len,len,2);
Pauli_Y_array = zeros(vec_len,vec_len,len,2);

for i = 1:len
    %% 1 for positive
    Pauli_X_array(:,:,i,1) = 0.5 * (tensor_exp(bit_I,len) + kron(kron(tensor_exp(bit_I,i-1),bit_X),tensor_exp(bit_I,len-i)));
    Pauli_Z_array(:,:,i,1) = 0.5 * (tensor_exp(bit_I,len) + kron(kron(tensor_exp(bit_I,i-1),bit_Z),tensor_exp(bit_I,len-i)));
    Pauli_Y_array(:,:,i,1) = 0.5 * (tensor_exp(bit_I,len) + kron(kron(tensor_exp(bit_I,i-1),bit_Y),tensor_exp(bit_I,len-i)));
     %% 2 for negative
    Pauli_X_array(:,:,i,2) = 0.5 * (tensor_exp(bit_I,len) - kron(kron(tensor_exp(bit_I,i-1),bit_X),tensor_exp(bit_I,len-i)));
    Pauli_Z_array(:,:,i,2) = 0.5 * (tensor_exp(bit_I,len) - kron(kron(tensor_exp(bit_I,i-1),bit_Z),tensor_exp(bit_I,len-i)));
    Pauli_Y_array(:,:,i,2) = 0.5 * (tensor_exp(bit_I,len) - kron(kron(tensor_exp(bit_I,i-1),bit_Y),tensor_exp(bit_I,len-i)));
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

%% test pauli projectors
for i = 1:10 %projector only test
    projector_choice = randi(6,1,1);
    if projector_choice == 1  %+X
        bit_choice = randi(len,1,1);
        projector = Pauli_X_array(:,:,bit_choice,1);
        s=s.CH_pauli_proj('+X',bit_choice);
    elseif projector_choice == 2  %-X
        bit_choice = randi(len,1,1);
        projector = Pauli_X_array(:,:,bit_choice,2);
        s=s.CH_pauli_proj('-X',bit_choice);
    elseif projector_choice == 3  %+Z
        bit_choice = randi(len,1,1);
        projector = Pauli_Z_array(:,:,bit_choice,1);
        s=s.CH_pauli_proj('+Z',bit_choice);
    elseif projector_choice == 4  %-Z
        bit_choice = randi(len,1,1);
        projector = Pauli_Z_array(:,:,bit_choice,2);
        s=s.CH_pauli_proj('-Z',bit_choice);
    elseif projector_choice == 5  %+Y
        bit_choice = randi(len,1,1);
        projector = Pauli_Y_array(:,:,bit_choice,1);
        s=s.CH_pauli_proj('+Y',bit_choice);
    elseif projector_choice == 6  %-Y
        bit_choice = randi(len,1,1);
        projector = Pauli_Y_array(:,:,bit_choice,2);
        s=s.CH_pauli_proj('-Y',bit_choice);
    else
        fprintf('error invalid projector choice.\n');
    end
    state_vector = projector * state_vector;
    s_state_vec = CH2basis(s);
    fprintf('%dth projector %d test!\n',i,projector_choice);
    %disp(bit_choice);
    %disp(state_vector); disp(s_state_vec);
    s.pp_CH('ch');
    assert(approx_equal(state_vector,s_state_vec,0.000000001)); %% may need +- to account for rounding error...
    fprintf('%dth projector %d passed!\n',i,projector_choice);
end
