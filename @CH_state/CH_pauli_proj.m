%%
%%  stab_new = stab.CH_pauli_proj(g,q)
%%               
function stab_new = CH_pauli_proj(obj,neg,x_bit,z_bit)

    stab_new = CH_state(obj.len);
    stab_new.deepcopy(obj);

    t = obj.s;
    u = obj.s; 
    a = const.init_uint;
    d = mod(bitsum(bitand(x_bit,z_bit)),4); % imaginary part from Y's
    if neg
        d = bitset(d,2,~bitget(d,2));  % sign from choice of pauli
    end
    % Z part
    new_zz_bit = const.init_uint;
    for j = 1:obj.len
        if bitget(z_bit,j)
            new_zz_bit = bitxor(new_zz_bit,obj.G(j,1));
        end
    end


    % X part  
    new_x_bit = const.init_uint;
    new_xz_bit = const.init_uint;

    for j = 1:obj.len
        if bitget(x_bit,j)
            d = mod(d + obj.get_g(j),4);
            if bitget(parity(bitand(obj.M(j,1),new_x_bit)),1)
                d = bitset(d,2,~bitget(d,2));
            end
            new_x_bit = bitxor(new_x_bit,obj.F(j,1));
            new_xz_bit = bitxor(new_xz_bit,obj.M(j,1));
        end
    end

    x_bit = new_x_bit;
    z_bit = bitxor(new_zz_bit,new_xz_bit);

    % H part
    %second layer
    ZqL = bitand(x_bit,stab_new.v);
    XqR = bitand(z_bit,stab_new.v);
    % first layer
    ZqR = bitxor(XqR,z_bit);
    XqL = bitxor(ZqL,x_bit);

    neg = parity(bitand(ZqL,XqR));
    neg = parity(bitand(ZqR,t)) + neg;
    neg = parity(bitand(ZqL,t)) + neg;
    neg = bitget(neg,1);
    u = bitxor(XqR,t); 
    u = bitxor(XqL,u); 
    if neg
        d = bitset(d,2,~bitget(d,2));
    end

    if t ~= u
        % todo: refactor superpos2circuit to take 2.^(0.5) out 
        superpos2circuit(t,u,a,d,stab_new);
        stab_new.w = 0.5 * stab_new.w * 2.^(0.5);
    else
        stab_new.w = stab_new.w * 0.5 * (1+(1i).^double(d));
    end
    if stab_new.w ~= 0
        stab_new.w = stab_new.w/norm(stab_new.w);
    end
 
end
