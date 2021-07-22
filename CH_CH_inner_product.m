%%
%%  inner_product = CH_CH_inner_product(astab,bstab)                                                        
%%  computes inner product <astab|bstab>, where |astab>,|bstab> are stabilizer states of same bit length
%%
function inner_product = CH_CH_inner_product(astab,bstab)
    % deepcopy of bstab, so we don't mutate astab or bstab
    abstab = CH_state(bstab.len);
    abstab.deepcopy(bstab);
    
    % compute Ucab = (Uca^*)Ucb tableau 
    % Z part
    for i = 1:abstab.len
        Gab_i = const.init_uint;
        GTa_i = astab.GT(i,1);
        for j = 1:abstab.len
            if bitget(GTa_i,j)
                Gab_i = bitxor(Gab_i,bstab.G(j,1));
            end
        end
        abstab.G(i,1) = Gab_i;
    end
    
    % X part  
    for i = 1:abstab.len
        Fab_i = const.init_uint;
        Mab_i = const.init_uint;
        mu_i = astab.get_gT(i); % from product of mu_i's
        mu_i_sign = const.init_uint; % from pulling out Z gates
        FTa_i = astab.FT(i,1);
        MTa_i = astab.MT(i,1);
        for j = 1:abstab.len
            if bitget(FTa_i,j)
                mu_i_sign = bitxor(mu_i_sign,parity(bitand(bstab.M(j,1),Fab_i)));
                Fab_i = bitxor(Fab_i,bstab.F(j,1));
                Mab_i = bitxor(Mab_i,bstab.M(j,1));
                mu_i = mod(mu_i+bstab.get_g(j),4);
            end
            if bitget(MTa_i,j)
                Mab_i = bitxor(Mab_i,bstab.G(j,1));
            end
        end

        abstab.F(i,1) = Fab_i;
        abstab.M(i,1) = Mab_i;
        if mu_i_sign
            mu_i = mod(mu_i+2,4);
        end
        abstab.set_g(i,mu_i);
    end
    abstab.w = abstab.w * conj(astab.w);

    % apply Uha 
    for i = 1:abstab.len
        if bitget(astab.v,i)
            abstab.CH_gate('HL',i);
        end
    end

    inner_product = CH_basis_inner_product(astab.s,abstab);
end
