%% 
%% generate x string of commuter pauli for set of generators gen_array
%% assumes gen(i) in {I,Z}^len for all i
%% 
function x_bits = get_commuter(gen_array,leading_bits,len) 
    bit_max = 2.^len;
    x_bits = cast(randi(bit_max,1,1) - 1,const.typecast_str);
    x_bits = bitxor(x_bits,leading_bits);
    x_bits = bitand(x_bits,bitcmp(leading_bits));
    for i = 1:len
        if gen_array(i) == 0
            continue;
        else
            bit_i = bitsum(bitset(bitand(gen_array(i),x_bits),i,0));
            x_bits = bitset(x_bits,i,bit_i);
        end
    end
end
