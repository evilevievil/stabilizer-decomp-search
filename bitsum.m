%% calculates the bitsum (weight) of x in binary
function bit_sum = bitsum(x)
    bit_sum = const.init_uint;
    for i = 1:const.init_max_qubits
        bit_sum = bit_sum + bitget(x,i);
    end
end
