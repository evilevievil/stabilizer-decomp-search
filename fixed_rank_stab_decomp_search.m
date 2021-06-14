%% stab_decomp_search(a)
%% main algorithm that searches for low-rank stabilizer decomposition of arbitrary state a

function stab_decomp = fixed_rank_stab_decomp_search(a,len,decomp_len,b_init,b_final,sa_max_step,walk_max_step)
    %%%%%%%%%%% init %%%%%%%%%%%
    vec_len = 2.^len;
    b = b_init; % beta
    step_ratio = (b_final / b_init).^(1/sa_max_step);
    obj_val = 0;
    
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

    for i = 1:sa_max_step
        result = random_walk(reverse_formatted_a,len,stab_decomp,decomp_len,walk_max_step,b,obj_val);
        [obj_val,new_stab_decomp] = result;
        b = b * step_ratio;
        if abs(obj_val-1) < 0.0000000000000001
            fprintf('found decomp with rank %d\n',decomp_len);
            break;
        end
    end
end

%% returns [projective norm, update stab decomp]
function result = random_walk(a,len,stab_decomp,decomp_len,walk_max_step,b,prev_walk_result)
    %% todo: try and optimize pauli proj (we can reduce # of failed pauli projections)
    prev_walk_result = 0;
    for i = 1:walk_max_step
        result = pauli_update(a,len,stab_decomp,decomp_len);
        curr_walk_result = result(1);
        new_stab_decomp = result(2)
        if abs(curr_walk_result-1) < 0.0000000000000001
            stab_decomp = new_stab_decomp; %% to do optimize by avoiding rededundant work
            prev_walk_result = curr_walk_result;
            break; %% found decomp
        elseif curr_walk_result > prev_walk_result
            stab_decomp = new_stab_decomp;
            prev_walk_result = curr_walk_result;
        else
            accept_pr = exp(-b*(prev_walk_result-curr_walk_result));
            if accept_pr <= rand()
                stab_decomp = new_stab_decomp;
                prev_walk_result = curr_walk_result;
            end
        end
    end
    result = [curr_walk_result,stab_decomp];
end

%% 
function result = pauli_update(reverse_formatted_a,len,stab_decomp,decomp_len)
    new_obj_val = -1;
    new_stab_decomp = stab_decomp;
    while new_obj_val == -1
        state_choice = randi(decomp_len,1,1);
        new_stab_decomp(state_choice) = stab_decomp(state_choice).CH_pauli_proj('rand',-1);
        new_obj_val = CH_decomp_project(reverse_formatted_a,new_stab_decomp,len,decomp_len);
    end
    result = [new_obj_val,new_stab_decomp];
end
