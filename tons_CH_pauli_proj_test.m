%%
%%
%% globals
len = 8;
vec_len = 2.^len;
bit_X = [0,1;1,0];
bit_Z = [1,0;0,-1];
bit_Y = [0,-1i;1i,0];
bit_I = [1,0;0,1];
big_I = tensor_exp(bit_I,len);



for k=1:100

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
        projector = 1;
        bit_choice = randi(4,len,1);
        sign_choice = randi(2,1,1);
        sign_choice = sign_choice-1;
        x_bits = const.init_uint;
        z_bits = const.init_uint;
        for j = 1:len
            if bit_choice(j,1) == 1 % I  
                projector = kron(projector,bit_I);
            elseif bit_choice(j,1) == 2 % X
                x_bits = bitset(x_bits,j,1);
                projector = kron(projector,bit_X);
            elseif bit_choice(j,1) == 3 % Z
                z_bits = bitset(z_bits,j,1);
                projector = kron(projector,bit_Z);
            else % Y
                x_bits = bitset(x_bits,j,1);
                z_bits = bitset(z_bits,j,1);
                projector = kron(projector,bit_Y);
            end
        end
        if sign_choice
            projector = 0.5 * (big_I - projector);
        else
            projector = 0.5 * (big_I + projector);
        end
        s = s.CH_pauli_proj(sign_choice,x_bits,z_bits);
    
        state_vector = projector * state_vector;
        if dot(state_vector,state_vector) ~= 0
            state_vector = state_vector/(dot(state_vector,state_vector)).^0.5;
        end
        s_state_vec = CH2basis(s);
        fprintf('%dth projector test!\n',i);
        disp(bit_choice);
        disp(sign_choice);
        disp(state_vector); 
        disp(s_state_vec);
        s.pp_CH('ch');
        assert(approx_equal(state_vector,s_state_vec,0.000000001)); %% may need +- to account for rounding error...
        fprintf('%dth projector passed!\n',i);
    end

end
