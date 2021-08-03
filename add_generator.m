%% 
%% either init or add new generator to generator array gen_array
%% assumes gen(i) in {I,Z}^len for all i
%% 
function [new_gen_array,new_leading_bits,new_gen_bits] = add_generator(gen_array,leading_bits,len) 
    new_leading_bit = 0;
    new_gen_array = gen_array;
    while true
        new_leading_bit = randi(len,1,1);
        if new_gen_array(new_leading_bit) > 0
            continue;
        else
            bit_max = 2.^new_leading_bit;
            new_gen_bits = bitset(cast(randi([2,bit_max],1,1) - 1,const.typecast_str),new_leading_bit,1);
            new_gen_array(new_leading_bit) = new_gen_bits;
            new_leading_bits = bitset(leading_bits,new_leading_bit,1);
            break;
        end
    end
end
