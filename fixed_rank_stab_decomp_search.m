%% stab_decomp_search(a)
%% main algorithm that searches for low-rank stabilizer decomposition of arbitrary state a

function stab_decomp = fixed_rank_stab_decomp_search(a,len,decomp_len,b_init,b_final,sa_max_step,walk_max_step)
    %%%%%%%%%%% init %%%%%%%%%%%
    vec_len = 2.^len;
    b = b_init; % beta
    step_ratio = (b_final - b_init)/sa_max_step;
    %step_ratio = (b_final / b_init).^(1/sa_max_step);
    
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
        stab_decomp(i).CH_init('zero');
        stab_decomp(i).CH_init('rand');
    end
    %load('./stab_decomp.mat','stab_decomp');
    %stab_decomp(1).CH_gate('HL',1);
    %stab_decomp(1).CH_gate('CXL',[1,2]);
    %stab_decomp(1).CH_gate('SL',1);
    %stab_decomp(2).s = 1;
    %stab_decomp(2).CH_gate('HL',2);
    %stab_decomp(2).CH_gate('CXL',[2,1]);
    
    obj_val = CH_decomp_project(reverse_formatted_a,stab_decomp,len,decomp_len);
    fprintf('print initial states\n');
    % print initial states
    %disp(obj_val);
    %for j = 1:decomp_len
    %    fprintf('decomp state %d\n',j);
    %    stab_decomp(j).pp_CH('basis');
    %    %stab_decomp(j).pp_CH('ch');
    %end
    
    for i = 1:sa_max_step
        [new_obj_val,new_stab_decomp] = random_walk(reverse_formatted_a,len,stab_decomp,decomp_len,walk_max_step,b,obj_val);
        
        if abs(new_obj_val-1) < 0.000000001
            stab_decomp = new_stab_decomp; %% to do optimize by avoiding rededundant work
            obj_val = new_obj_val;
            disp(obj_val);
            fprintf('FOUND!!!!\n');
            break; %% found decomp
        %elseif abs(obj_val-1) < 0.000000001
            disp(obj_val);
        %    fprintf('FOUND!!!!\n');
        %    break; %% found decomp
        elseif new_obj_val > obj_val
            stab_decomp = new_stab_decomp;
            obj_val = new_obj_val;
            disp(obj_val);
            %for j = 1:decomp_len
            %    fprintf('decomp state %d\n',j);
            %    stab_decomp(j).pp_CH('basis');
            %    %stab_decomp(j).pp_CH('ch');
            %end
        else
            accept_pr = exp(-b*(obj_val-new_obj_val));
            %fprintf('accept pr:%d,b:%d\n',accept_pr,b);
            if accept_pr >= rand() 
                stab_decomp = new_stab_decomp;
                obj_val = new_obj_val;
                disp(obj_val);
                %for j = 1:decomp_len
                %    fprintf('decomp state %d\n',j);
                %    stab_decomp(j).pp_CH('basis');
                %    %stab_decomp(j).pp_CH('ch');
                %end
            end
        end
        
        b = b + step_ratio;
        %disp(obj_val);
        %if obj_val > 0.90
            %disp(obj_val);
            %for j = 1:decomp_len
            %    fprintf('decomp state %d\n',j);
            %    stab_decomp(j).pp_CH('basis');
            %    %stab_decomp(j).pp_CH('ch');
            %end
        %end
    end
    disp(obj_val);
    for j = 1:decomp_len
        fprintf('decomp state %d\n',j);
        stab_decomp(j).pp_CH('basis');
        %stab_decomp(j).pp_CH('ch');
    end
end

%% returns [projective norm, update stab decomp]
function [x,y] = random_walk(a,len,stab_decomp,decomp_len,walk_max_step,b,obj_val)
    %% todo: try and optimize pauli proj (we can reduce # of failed pauli projections)
    for i = 1:walk_max_step
       % if obj_val < 1
            [curr_walk_result,stab_decomp] = pauli_update(a,len,stab_decomp,decomp_len);
       % else
       %     [curr_walk_result,stab_decomp] = Uc_update(a,len,stab_decomp,decomp_len);
       % end
        %fprintf('rand walk end: %d\n',i);
    end
    x = curr_walk_result;
    y = stab_decomp;
end

%% 
function [x,y] = pauli_update(reverse_formatted_a,len,stab_decomp,decomp_len)
    new_obj_val = -1;
    new_stab_decomp = stab_decomp;
    while new_obj_val == -1
        state_choice = randi(decomp_len,1,1);
        [sign_choice,x_bits,z_bits] = random_pauli_bits(len,len); 
        new_stab_decomp(state_choice) = stab_decomp(state_choice).CH_pauli_proj(sign_choice,x_bits,z_bits);
        new_obj_val = CH_decomp_project(reverse_formatted_a,new_stab_decomp,len,decomp_len);
    end
    %fprintf('pauli project success!\n');
    x = new_obj_val;
    y = new_stab_decomp;
    %for i =1:decomp_len
    %   disp(CH_CH_inner_product(new_stab_decomp(i),new_stab_decomp(i)));
    %end
end

function [x,y] = Uc_update(reverse_formatted_a,len,stab_decomp,decomp_len)
    new_obj_val = -1;
    new_stab_decomp = stab_decomp;
    while new_obj_val == -1
        state_choice = randi(decomp_len,1,1);
        new_stab_decomp(state_choice) = stab_decomp(state_choice).CH_gate('rand',-1);
        new_obj_val = CH_decomp_project(reverse_formatted_a,new_stab_decomp,len,decomp_len);
    end
    %fprintf('pauli project success!\n');
    x = new_obj_val;
    y = new_stab_decomp;
    %for i =1:decomp_len
    %   disp(CH_CH_inner_product(new_stab_decomp(i),new_stab_decomp(i)));
    %end
end

%% 
function [sign_choice,x_bits,z_bits] = random_pauli_bits(len,proj_num)
    bit_choice = randi(4,len,1);
    proj_choice = randperm(len,proj_num);
    sign_choice = randi(2,1,1);
    sign_choice = sign_choice-1;
    x_bits = uint8(0);
    z_bits = uint8(0);
    for j = 1:proj_num
        bit_loc = proj_choice(j);
        if bit_choice(bit_loc,1) == 1 % I  
            continue;
        elseif bit_choice(bit_loc,1) == 2 % X
            x_bits = bitset(x_bits,bit_loc,1);
        elseif bit_choice(bit_loc,1) == 3 % Z
            z_bits = bitset(z_bits,bit_loc,1);
        else % Y
            x_bits = bitset(x_bits,bit_loc,1);
            z_bits = bitset(z_bits,bit_loc,1);
        end
    end
end
