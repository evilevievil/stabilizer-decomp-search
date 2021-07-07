%% stab_decomp_search(a)
%% main algorithm that searches for low-rank stabilizer decomposition of arbitrary state a

%% main search function
function stab_decomp = fixed_rank_stab_decomp_search(a,len,decomp_len,b_init,b_final,sa_max_step,walk_max_step)
    %%%%%%%%%%% init %%%%%%%%%%%
    %rng(89);
    %rng(6); % t_5_5 h_6_7
    rng(876); % h_6_7 pt2
    b = b_init; 
    step_ratio = (b_final - b_init)/sa_max_step;
    reverse_formatted_a = reverse_format_amp(a,len);
    new_obj_val = 0;

    %% init stab_decomp
    prev_data = load('H_6_7_0.9539.mat'); % load from saved data
    for i= 1:decomp_len
        stab_decomp(i) = CH_state(len);
        stab_decomp(i).CH_init('zero');
        %stab_decomp(i).CH_init('rand');
        stab_decomp(i).deepcopy(prev_data.ans(i)); % load from saved data
    end
    
    [obj_val,G,a_stab_array] = CH_decomp_project(reverse_formatted_a,stab_decomp,len,decomp_len);
    %% BUG??? check initial states are linearly independent 
    assert(obj_val~=-1);

    %%%%%%%%%%% search %%%%%%%%%%%
    fprintf('Search Start!\n');
    disp(obj_val);
    
    % SA cooling loop
    for i = 1:sa_max_step

        for j = 1:walk_max_step
            [new_obj_val,new_stab_decomp,new_G,new_a_stab_array] = pauli_update(reverse_formatted_a,G,a_stab_array,len,stab_decomp,decomp_len);
            if abs(new_obj_val-1) < 0.000000001
                stab_decomp = new_stab_decomp; 
                obj_val = new_obj_val;
                G = new_G;
                a_stab_array = new_a_stab_array;
                disp(obj_val);
                fprintf('FOUND!!!!\n');
                break; 
            elseif new_obj_val > obj_val
                stab_decomp = new_stab_decomp; 
                obj_val = new_obj_val;
                G = new_G;
                a_stab_array = new_a_stab_array;
                disp(new_obj_val);
            else
                accept_pr = exp(-b*(obj_val-new_obj_val));
                if accept_pr >= rand() 
                    stab_decomp = new_stab_decomp; 
                    obj_val = new_obj_val;
                    G = new_G;
                    a_stab_array = new_a_stab_array;
                    disp(new_obj_val);
                end
                %disp(new_obj_val);
            end
        end
        %disp(obj_val);
        
        if abs(new_obj_val-1) < 0.000000001
            break;
        end
        b = b + step_ratio;
    end
    disp(obj_val);
    for j = 1:decomp_len
        fprintf('decomp state %d\n',j);
        stab_decomp(j).pp_CH('basis');
    end
end

%% 
function [new_obj_val,new_stab_decomp,new_G,new_a_stab_array] = pauli_update(reverse_formatted_a,G,a_stab_array,len,stab_decomp,decomp_len)
    state_choice = 1;
    new_obj_val = -1;
    new_stab_decomp = stab_decomp;
    while new_obj_val == -1 || (CH_CH_inner_product(new_stab_decomp(state_choice),stab_decomp(state_choice))==1)
    %while new_obj_val == -1 
        new_stab_decomp(state_choice) = stab_decomp(state_choice);
        state_choice = randi(decomp_len,1,1);
        [sign_choice,x_bits,z_bits] = random_pauli_bits(len); 
        new_stab_decomp(state_choice) = stab_decomp(state_choice).CH_pauli_proj(sign_choice,x_bits,z_bits);
        [new_obj_val,new_G,new_a_stab_array] = CH_decomp_project_memoize(reverse_formatted_a,new_stab_decomp,G,a_stab_array,state_choice,len,decomp_len);
    end
end

%% 
% todo: generate random number as bit string to optimize
function [sign_choice,x_bits,z_bits] = random_pauli_bits(len)
    bit_choice = randi(3,len,1);
    sign_choice = randi(2,1,1);
    sign_choice = sign_choice-1;
    x_bits = const.init_uint;
    z_bits = const.init_uint;
    for j = 1:len
        if bit_choice(j,1) == 1 % I  
            continue;
        elseif bit_choice(j,1) == 2 % X
            x_bits = bitset(x_bits,j,1);
        elseif bit_choice(j,1) == 3 % Z
            z_bits = bitset(z_bits,j,1);
        else % Y
            %x_bits = bitset(x_bits,j,1);
            %z_bits = bitset(z_bits,j,1);
        end
    end
end
