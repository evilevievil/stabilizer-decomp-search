%% magic_state_vec computes state vector array for len qubit magic_state given type 
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
    case 'r_1_3'
        code_words = cast([bin2dec('00000000'),bin2dec('00001111'),bin2dec('01010101'),bin2dec('01011010'),...
                           bin2dec('10101010'),bin2dec('10100101'),bin2dec('11111111'),bin2dec('11110000'),...
                           bin2dec('00110011'),bin2dec('00111100'),bin2dec('01100110'),bin2dec('01101001'),...
                           bin2dec('10011001'),bin2dec('10010110'),bin2dec('11001100'),bin2dec('11000011')],const.typecast_str);
        magic_state = 0;
        for i=1:16
            magic_state = magic_state + bin2codestate(code_words(i),8);
        end
        magic_state = magic_state/norm(magic_state);
    end
end

function state_vec = bin2codestate(b,len)
    bit_T = [2.^-0.5; 0.5*(1+1i)];
    bit_Tt = [2.^-0.5; -0.5*(1+1i)];
    state_vec = 1;
    for i=1:len
        if bitget(b,i)
            state_vec = kron(state_vec,bit_Tt);
        else
            state_vec = kron(state_vec,bit_T);
        end      
    end
end