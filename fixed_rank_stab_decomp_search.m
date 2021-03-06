%% stab_decomp_search(a)
%% main algorithm that searches for low-rank stabilizer decomposition of target state a

%% main search function
function stab_decomp = fixed_rank_stab_decomp_search(a,len,decomp_len,b_init,b_final,sa_max_step,walk_max_step)
    %%%%%%%%%%% log graph data %%%%%%%%%%%
    %search_status = 0;
    %all_x = 1:(50*1000);
    %all_y = zeros(1,50*1000);
    %SA_x = 1:50;
    %SA_y = zeros(1,50);

    %%%%%%%%%%% init %%%%%%%%%%%
    rng(1); 
    b = b_init; 
    step_ratio = (b_final - b_init)/sa_max_step;
    reverse_formatted_a = reverse_format_amp(a,len);
    new_obj_val = 0;

    %% init stab_decomp
    prev_data = load('12_gen_cat_0.51.mat'); % load from saved data
    for i= 1:decomp_len
        stab_decomp(i) = CH_state(len);
        stab_decomp(i).CH_init('zero');
        %stab_decomp(i).CH_init('rand');
        stab_decomp(i).deepcopy(prev_data.ans(i)); % load from saved data
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
        %for j = 1:walk_steps(i)
        for j = 1:walk_max_step
            [new_obj_val,new_stab_decomp,new_G,new_a_stab_array] = pauli_update(reverse_formatted_a,G,a_stab_array,len,stab_decomp,decomp_len);
            if abs(new_obj_val-1) < 0.000000001 % account for rounding error
                stab_decomp = new_stab_decomp; 
                obj_val = new_obj_val;
                G = new_G;
                a_stab_array = new_a_stab_array;
                disp(obj_val);
                % search_status = 1; for plot data
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
            end
            %all_y(1,walk_max_step*(i-1)+j) = obj_val; % log graph data
        end

        if abs(new_obj_val-1) < 0.000000001
            break;
        end
        %SA_y(1,i) = obj_val; % log graph data
        b = b + step_ratio;
    end

    disp(obj_val);
    for j = 1:decomp_len
        fprintf('decomp state %d\n',j);
        stab_decomp(j).pp_CH('basis');
    end
    %save('all_x.mat','all_x'); % log graph data
    %save('all_y.mat','all_y'); % log graph data
    %save('SA_x.mat','SA_x'); % log graph data
    %save('SA_y.mat','SA_y'); % log graph data
end

%% perform 1 walk step by updating stab decomp array
function [new_obj_val,new_stab_decomp,new_G,new_a_stab_array] = pauli_update(reverse_formatted_a,G,a_stab_array,len,stab_decomp,decomp_len)
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
