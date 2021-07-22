function par = parity(x)
    par = const.init_uint;
    for i = 1:const.init_max_qubits
        par = bitxor(par,bitget(x,i));
    end
end
