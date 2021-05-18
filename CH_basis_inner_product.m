%%
%%  inner_product = CH_basis_inner_product(a,stab)                                                        
%%  computes inner product <stab_i|stab_j>, where both are stabilizer states of bit length len
%%
function inner_product = CH_basis_inner_product(x,stab)
    t = uint8(0);
    u = uint8(0);
    mu = uint8(0);
    mu_sign = uint8(0);
    % compute t for Z, u for X, mu for phase
    for i = 1:stab.len
        if bitget(x,i)
            t = bitxor(t,stab.M(i,1));
            u = bitxor(u,stab.F(i,1));
            mu = mu + stab.get_g(i);
        end
    end
    mu = mod(mu,4); 
    % compute additional -1 mu_sign from swapping all Z's to the left 
    %for i = 1:stab.len
    %    if bitget(t,i)
    %        for j = 1 : i
    %            mu_sign = bitxor(mu_sign,bitget(u,j));
    %        end
    %    end
    %end
    mu_sign = parity(bitand(t,u));
    % compute neg and sign from applying Uc*XUcUh on s
    s = stab.s;
    v = stab.v;
    v_len = bitsum(v);
    neg = bitxor(parity(bitand(v,bitand(u,s))),mu_sign);
    one = bitand(bitcmp(v),bitxor(bitcmp(u),s));
    if one == bitcmp(v)
        inner_product = stab.w * (-1).^(double(neg)) * (1i).^(double(mu)) * 2.^(-0.5 * double(v_len));
    else
        inner_product = double(0);
    end
end