%%
%%  inner_product = CH_basis_inner_product(x,stab)                                                        
%%  computes inner product <x_str|stab> where x_str is the stab.len bit binary representation of x
%%
function inner_product = CH_basis_inner_product(x,stab)
    % init zero bit strings
    t = const.init_uint;
    u = const.init_uint;
    mu = const.init_uint;
    mu_sign = const.init_uint;

    % compute t for Z, u for X, mu for phase
    for i = 1:stab.len
        if bitget(x,i)
            t = bitxor(t,stab.M(i,1));
            u = bitxor(u,stab.F(i,1));
            mu = mu + stab.get_g(i);
            mu_sign = mu_sign + parity(bitand(u,stab.M(i,1)));
        end
    end
    mu = mod(mu,4); 

    % compute neg and sign from applying Uc*XUcUh on s
    s = stab.s;
    v = stab.v;
    v_len = bitsum(v);
    neg = parity(bitand(v,bitand(u,s)))+mu_sign;
    one = bitand(bitcmp(v),bitxor(bitcmp(u),s));

    if one == bitcmp(v)
        inner_product = stab.w * (-1).^(double(neg)) * (1i).^(double(mu)) * 2.^(-0.5 * double(v_len));
    else
        inner_product = double(0);
    end
end
