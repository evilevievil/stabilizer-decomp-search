%%
%%  inner_product = CH_CH_inner_product(astab,bstab)                                                        
%%  computes inner product <astab|bstab>, where |astab>,|bstab> are stabilizer states of same bit length
%%
% computing sign bit is slow..
function inner_product = CH_CH_inner_product(astab,bstab)
    abstab = CH_state(bstab.len);
    abstab.deepcopy(bstab);
    
    % compute Ucab = (Uca^*)Ucb tableau 
    % Z part
    for i = 1:abstab.len
        Gab_i = uint8(0);
        for j = 1:abstab.len
            if astab.get_G(i,j)
                Gab_i = bitxor(Gab_i,bstab.G(j,1));
            end
        end
        abstab.G(i,1) = Gab_i;
    end
    
    % X part  (problematic...)
    for i = 1:abstab.len
        Fab_i = uint8(0);
        Mab_i = uint8(0);
        mu_i = uint8(0);
        mu_i_sign = uint8(0);
        for j = 1:abstab.len
            if astab.get_F(i,j)
                Fab_i = bitxor(Fab_i,bstab.F(j,1));
            end
            if astab.get_M(i,j)
                Mab_i = bitxor(Mab_i,bstab.M(j,1));
            end
        end
        abstab.F(i,1) = Fab_i;
        abstab.M(i,1) = Mab_i;
    end

    % apply Uha 
    for i = 1:abstab.len
        if bitget(astab.v,i)
            abstab.CH_gate('HL',i);
        end
    end

    inner_product = CH_basis_inner_product(astab.s,abstab);
end


