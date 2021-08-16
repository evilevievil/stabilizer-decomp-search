%% stab_decomp_search(a)
%% main algorithm that searches for low-rank stabilizer decomposition of target state a

%% main search function
function [stab_decomp,gen_array,k] = magic_state_decomp_search_v2(len,decomp_len,b_init,b_final,sa_max_step,walk_max_step,max_idle_steps,max_k,max_k_iter)

    %%%%%%%%%%% init %%%%%%%%%%%
    rng(265436); 
    best_so_far = 0;
    if max_k <0
       max_k = floor(0.5*(len - (4*log2(decomp_len)/log2(3)))); 
    end
    a = magic_state_vec('T',len);
    b = b_init; 
    step_ratio = (b_final - b_init)/sa_max_step;
    reverse_formatted_a = reverse_format_amp(a,len);
    new_obj_val = 0;
    idle_step = 0;
    gen_array = zeros(len,1);
    gen_array = cast(gen_array,const.typecast_str);
    leading_bits = cast(0,const.typecast_str);
    target_obj_val = 1;
    k = 0;
    k_iter = 0;

    %% init stab_decomp
    %prev_data = load('catT_8_5_0.9133.mat'); % load from saved data
    for i= 1:decomp_len
        stab_decomp(i) = CH_state(len);
        stab_decomp(i).CH_init('zero');
        stab_decomp(i).CH_init('rand');
        %stab_decomp(i).deepcopy(prev_data.stab_decomp(i)); % load from saved data
    end

    [obj_val,G,a_stab_array] = CH_decomp_project(reverse_formatted_a,stab_decomp,len,decomp_len);
    %% check initial states are linearly independent 
    assert(obj_val~=-1);
    
    %%%%%%%%%%% search %%%%%%%%%%%
    fprintf('Search Start!\n');
    disp(obj_val);
    
    % SA cooling loop
    for i = 1:sa_max_step
        fprintf('SA %dth step\n',i);
        idle_walk = 1;
        for j = 1:walk_max_step
            [new_obj_val,new_stab_decomp,new_G,new_a_stab_array] = pauli_update(reverse_formatted_a,G,a_stab_array,len,stab_decomp,decomp_len,gen_array,leading_bits);
            if abs(new_obj_val-target_obj_val) < 0.000000001 % account for rounding error
                stab_decomp = new_stab_decomp; 
                obj_val = new_obj_val;
                G = new_G;
                a_stab_array = new_a_stab_array;
                idle_walk = 0;
                disp(obj_val);
                fprintf('FOUND!!!!\n');
                break; 
            elseif new_obj_val > obj_val
                stab_decomp = new_stab_decomp; 
                obj_val = new_obj_val;
                G = new_G;
                a_stab_array = new_a_stab_array;
                if new_obj_val > best_so_far
                    best_so_far = new_obj_val;
                    idle_walk = 0;
                end
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
            end
        end

        if abs(new_obj_val-target_obj_val) < 0.000000001
            break;
        end
        
        %% heuristic for adding generators
        if idle_walk
            idle_step = idle_step + 1;
        else
            idle_step = 0;
        end

        if idle_step >= max_idle_steps
            idle_step = 0;
            if k+1 > max_k
                max_idle_steps = 100;
            elseif k_iter < max_k_iter && k>0
                [a,gen_array,leading_bits] = generator_search_v2(a,len,gen_array,leading_bits,stab_decomp,decomp_len,0,b,k);
                reverse_formatted_a = reverse_format_amp(a,len);
                [obj_val,G,a_stab_array] = CH_decomp_project(reverse_formatted_a,stab_decomp,len,decomp_len);
                k_iter = k_iter + 1;
                best_so_far = obj_val;
            else
                [a,gen_array,leading_bits] = generator_search_v2(a,len,gen_array,leading_bits,stab_decomp,decomp_len,1,b,k);
                reverse_formatted_a = reverse_format_amp(a,len);
                [obj_val,G,a_stab_array] = CH_decomp_project(reverse_formatted_a,stab_decomp,len,decomp_len);
                target_obj_val = target_obj_val/2;
                k = k + 1;
                k_iter = 0;
                best_so_far = obj_val;
            end
        end

        b = b + step_ratio;
    end

    disp(obj_val);
    for j = 1:decomp_len
        fprintf('decomp state %d\n',j);
        stab_decomp(j).pp_CH('basis');
    end
end

%% perform 1 walk step by updating stab decomp array
function [new_obj_val,new_stab_decomp,new_G,new_a_stab_array] = pauli_update(reverse_formatted_a,G,a_stab_array,len,stab_decomp,decomp_len,gen_array,leading_bits)
    state_choice = 1;
    new_obj_val = -1;
    new_stab_decomp = stab_decomp;
    while new_obj_val == -1 || (CH_CH_inner_product(new_stab_decomp(state_choice),stab_decomp(state_choice))==1)
        new_stab_decomp(state_choice) = stab_decomp(state_choice);
        state_choice = randi(decomp_len,1,1);
        [sign_choice,x_bits,z_bits] = random_pauli_bits(len); 
        new_stab_decomp(state_choice) = stab_decomp(state_choice).CH_pauli_proj(sign_choice,x_bits,z_bits);
        [new_obj_val,new_G,new_a_stab_array] = CH_decomp_project_memoize(reverse_formatted_a,new_stab_decomp,G,a_stab_array,state_choice,len,decomp_len);
    end
end

%% random pauli projector selector
function [sign_choice,x_bits,z_bits] = random_pauli_bits(len)
    bit_max = 2.^len;
    sign_choice = randi(2,1,1);
    sign_choice = sign_choice-1;
    x_bits = cast(randi(bit_max,1,1) - 1,const.typecast_str);
    z_bits = cast(randi(bit_max,1,1) - 1,const.typecast_str);
end
