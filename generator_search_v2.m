%%
%% generator_search: searches for magic code generators given current stab decomp and inver temperature
%%
function [a,gen_array,leading_bits] = generator_search_v2(a,len,gen_array,leading_bits,stab_decomp,decomp_len,isupdate,b,k)
    obj_val = -1;
    new_gen_array = gen_array;
    new_leading_bits = leading_bits;
    new_a = a;
    old_a = a;
    old_gen_array = gen_array;
    old_leading_bits = leading_bits;
    fprintf('generator search start\n');
    if isupdate
        for step = 1:1000
            [new_a,new_gen_array,new_leading_bits] = gen_update(old_a,len,old_gen_array,old_leading_bits);
            reverse_formatted_a = reverse_format_amp(new_a,len);
            [new_obj_val,~,~] = CH_decomp_project(reverse_formatted_a,stab_decomp,len,decomp_len);
            if new_obj_val > obj_val
                gen_array = new_gen_array;
                leading_bits = new_leading_bits;
                obj_val = new_obj_val;
                a = new_a;
            end 
        end
    else
        for step = 1:1000
            [new_a,new_gen_array,new_leading_bits] = gen_walk(old_a,len,old_gen_array,old_leading_bits,k);
            reverse_formatted_a = reverse_format_amp(new_a,len);
            [new_obj_val,~,~] = CH_decomp_project(reverse_formatted_a,stab_decomp,len,decomp_len);
            accept_pr = exp(-b*(obj_val-new_obj_val));
            if new_obj_val > obj_val || accept_pr >= rand() 
                gen_array = new_gen_array;
                leading_bits = new_leading_bits;
                obj_val = new_obj_val;
                a = new_a;
            end 
        end
    end
    fprintf('generator search end\n');
end

%% update generator array, objective value and project stab decomp
function [a,gen_array,leading_bits] = gen_walk(a,len,gen_array,leading_bits,k)
    a = magic_state_vec('T',len);
    zero_index = randi(k,1,1);
    bit_Z = [1,0;0,-1];
    bit_I = [1,0;0,1];
    len_I = tensor_exp(bit_I,len);
    curr_index = 1;
    for j = 1:len
        if bitget(leading_bits,j)
            if curr_index == zero_index
                leading_bits = bitset(leading_bits,j,0);
                gen_array(j) = 0;
            else
                projector_Z = 1;
                gen_bits = gen_array(j);
                for i=1:len
                    if bitget(gen_bits,i)
                        projector_Z = kron(projector_Z,bit_Z);
                    else
                        projector_Z = kron(projector_Z,bit_I);
                    end
                end
                a = 0.5 * (len_I + projector_Z) * a;
            end
            curr_index = curr_index + 1;
        end
    end
    [a,gen_array,leading_bits] = gen_update(a,len,gen_array,leading_bits);
end

%% update generator array, objective value and project stab decomp
function [a,gen_array,leading_bits] = gen_update(a,len,gen_array,leading_bits)
    new_gen_array = gen_array;
    new_leading_bits = leading_bits;
    bit_Z = [1,0;0,-1];
    bit_I = [1,0;0,1];
    len_I = tensor_exp(bit_I,len);
    new_a = zeros(len);
    while new_a' * new_a < 0.000000001
        projector_Z = 1;
        [new_gen_array,new_leading_bits,new_gen_bits] = add_generator(gen_array,leading_bits,len);
        for i=1:len
            if bitget(new_gen_bits,i)
                projector_Z = kron(projector_Z,bit_Z);
            else
                projector_Z = kron(projector_Z,bit_I);
            end
        end
        new_a = 0.5 * (len_I + projector_Z) * a;
    end
    a = new_a;
    gen_array = new_gen_array;
    leading_bits = new_leading_bits;
end

