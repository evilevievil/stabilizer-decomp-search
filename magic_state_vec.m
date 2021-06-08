function magic_state = magic_state_vec(type, len)
    bit_T = [2.^-0.5; 0.5*(1+1i)];
    switch type
    case 'T'
        magic_state = tensor_exp(bit_T,len);
    end
end