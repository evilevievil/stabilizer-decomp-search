%%%% get commuter test
%add_generator_test(6);
gen_update_test(6,3,2)

function get_commuter_test()
    gen_array = zeros(6,1);
    gen_array = cast(gen_array,const.typecast_str);
    gen_array(5) = cast(bin2dec('0000000000010001'),const.typecast_str);
    gen_array(3) = cast(bin2dec('0000000000000101'),const.typecast_str);
    leading_bits = cast(bin2dec('0000000000010100'),const.typecast_str);

    x_bits = get_commuter(gen_array,leading_bits,6);

    dec2bin(x_bits);

end

function add_generator_test(len)
    gen_array = zeros(len,1);
    gen_array = cast(gen_array,const.typecast_str);
    leading_bits = cast(0,const.typecast_str);

    for i = 1:len
        [gen_array,leading_bits,new_bit_str] = add_generator(gen_array,leading_bits,len);
        disp(dec2bin(new_bit_str));
    end

end

function gen_update_test(len,decomp_len,k)
    % empty generator array
    gen_array = zeros(len,1);
    gen_array = cast(gen_array,const.typecast_str);
    leading_bits = cast(0,const.typecast_str);
    % random stab decomp
    for i= 1:decomp_len
        stab_decomp(i) = CH_state(len);
        stab_decomp(i).CH_init('zero');
        stab_decomp(i).CH_init('rand');
    end
    a = magic_state_vec('T',len);
    reverse_formatted_a = reverse_format_amp(a,len);
    
    new_obj_val = -1;
    new_stab_decomp = stab_decomp;
    new_gen_array = gen_array;
    new_leading_bits = leading_bits;
    x_bits = cast(0,const.typecast_str);
    
    for i = 1:k
        new_obj_val = -1;
        while new_obj_val < 0
            disp(new_obj_val);
            [new_gen_array,new_leading_bits,new_gen_bits] = add_generator(gen_array,leading_bits,len);
            for i=1:decomp_len
                new_stab_decomp(i) = stab_decomp(i).CH_pauli_proj(0,x_bits,new_gen_bits);
            end
            [new_obj_val,~,~] = CH_decomp_project(reverse_formatted_a,new_stab_decomp,len,decomp_len);
        end
        gen_array = new_gen_array;
        leading_bits = new_leading_bits;
        stab_decomp = new_stab_decomp;
        obj_val = new_obj_val;
        disp(dec2bin(new_gen_bits));
    end
    %disp(dec2bin(gen_array));
end
