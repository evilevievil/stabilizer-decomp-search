%%
function [projection,G,a_stab_array] = CH_decomp_project_memoize(a,stab_array,G,a_stab_array,k,len,decomp_len)
    a_len = 2.^(len);
    projection = double(0);
    G_inv = zeros(decomp_len);

    % update G array
    for i = 1:decomp_len
        % todo: just use conjugate for lower triangular
        G(i,k) = CH_CH_inner_product(stab_array(i),stab_array(k));
        G(k,i) = CH_CH_inner_product(stab_array(k),stab_array(i));
    end


    % todo: can I optimize this?
    if (abs(G(k,k)-0) > 0.001) && (abs(det(G)-0) > 0.001)
        G_inv = inv(G);
        % update basis array
        a_stab_array(k) = 0;
        for j = 1:a_len
            a_stab_array(k) = a_stab_array(k)  + conj(a(j)) * CH_basis_inner_product(j-1,stab_array(k));
        end
        
        % calulate projection
        for i = 1:decomp_len
            for j = 1:decomp_len
                projection = projection + a_stab_array(i) * conj(a_stab_array(j)) * G_inv(i,j);
            end
        end
    else
        % Gram matrix is singular -> find a better stabilizer decomp
        projection = -1;
    end
end
