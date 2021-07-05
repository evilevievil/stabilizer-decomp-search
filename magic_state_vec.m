function magic_state = magic_state_vec(type, len)
    bit_T = [2.^-0.5; 0.5*(1+1i)];
    bit_H = [cos(pi/8); sin(pi/8)];
    bit_Tt = [2.^-0.5; -0.5*(1+1i)];
    switch type
    case 'T'
        magic_state = tensor_exp(bit_T,len);
    case 'H'
        magic_state = tensor_exp(bit_H,len);
    case 'catT'
        magic_state = 2.^(-0.5) * (tensor_exp(bit_T,len) + tensor_exp(bit_Tt,len));
    end
end