%% stab_decomp_search(a)
%% main algorithm that searches for low-rank stabilizer decomposition of target state a

%% main search function
function [stab_decomp,gen_array,k] = magic_code_decomp_search(len,decomp_len,b_init,b_final,sa_max_step,walk_max_step,max_idel_steps)

    %%%%%%%%%%% init %%%%%%%%%%%
    rng(0); 
    max_k = floor(0.5*(len - (4*log2(decomp_len)/log2(3))));
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
    disp(max_k);
    
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
            elseif abs(new_obj_val-target_obj_val) < abs(obj_val-target_obj_val) 
                stab_decomp = new_stab_decomp; 
                obj_val = new_obj_val;
                G = new_G;
                a_stab_array = new_a_stab_array;
                idle_walk = 0;
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

        if idle_step >= max_idel_steps
            if k+1 > max_k
                max_idel_steps = 100;
            else
                target_obj_val = target_obj_val/2;
                [obj_val,stab_decomp,G,a_stab_array,gen_array,leading_bits] = gen_update(reverse_formatted_a,len,stab_decomp,decomp_len,gen_array,leading_bits,target_obj_val);
                k = k + 1;
            end
        end

        b = b + step_ratio;
    end

    disp(obj_val);
    for j = 1:decomp_len
        fprintf('decomp state %d\n',j);
        stab_decomp(j).pp_CH('basis');
    end
    fprintf('k=%d\n',k);
    fprintf('max_k=%d\n',max_k);
end

%% update generator array, objective value and project stab decomp
function [obj_val,stab_decomp,G,a_stab_array,gen_array,leading_bits] = gen_update(reverse_formatted_a,len,stab_decomp,decomp_len,gen_array,leading_bits,target_obj_val)
    new_obj_val = -1;
    new_stab_decomp = stab_decomp;
    new_gen_array = gen_array;
    new_leading_bits = leading_bits;
    x_bits = cast(0,const.typecast_str);
    while new_obj_val < 0
        disp(new_obj_val);
        [new_gen_array,new_leading_bits,new_gen_bits] = add_generator(gen_array,leading_bits,len);
        for i=1:decomp_len
            new_stab_decomp(i) = stab_decomp(i).CH_pauli_proj(0,x_bits,new_gen_bits);
        end
        [new_obj_val,G,a_stab_array] = CH_decomp_project(reverse_formatted_a,new_stab_decomp,len,decomp_len);
    end
    gen_array = new_gen_array;
    leading_bits = new_leading_bits;
    stab_decomp = new_stab_decomp;
    obj_val = new_obj_val;
end

%% perform 1 walk step by updating stab decomp array
function [new_obj_val,new_stab_decomp,new_G,new_a_stab_array] = pauli_update(reverse_formatted_a,G,a_stab_array,len,stab_decomp,decomp_len,gen_array,leading_bits)
    state_choice = 1;
    new_obj_val = -1;
    new_stab_decomp = stab_decomp;
    while new_obj_val == -1 || (CH_CH_inner_product(new_stab_decomp(state_choice),stab_decomp(state_choice))==1)
        new_stab_decomp(state_choice) = stab_decomp(state_choice);
        state_choice = randi(decomp_len,1,1);
        [sign_choice,x_bits,z_bits] = random_pauli_bits(gen_array,leading_bits,len); 
        new_stab_decomp(state_choice) = stab_decomp(state_choice).CH_pauli_proj(sign_choice,x_bits,z_bits);
        [new_obj_val,new_G,new_a_stab_array] = CH_decomp_project_memoize(reverse_formatted_a,new_stab_decomp,G,a_stab_array,state_choice,len,decomp_len);
    end
end

%% random pauli projector selector
function [sign_choice,x_bits,z_bits] = random_pauli_bits(gen_array,leading_bits,len)
    bit_max = 2.^len;
    sign_choice = randi(2,1,1);
    sign_choice = sign_choice-1;
    x_bits = get_commuter(gen_array,leading_bits,len); 
    z_bits = cast(randi(bit_max,1,1) - 1,const.typecast_str);
end
