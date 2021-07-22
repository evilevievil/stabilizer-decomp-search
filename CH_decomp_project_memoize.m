%% [projection_value,G_new,a_stab_array_new] = 
%%        CH_decomp_project_memoize(target_state_vector, stab_array_new
%%                                  G_old,a_stab_array_old,update_index_k,len,decomp_len)
%% computes projection of target a onto stab_array given updated stab_array(k) and previous G and a_stab_array
function [projection,G,a_stab_array] = CH_decomp_project_memoize(a,stab_array,G,a_stab_array,k,len,decomp_len)
    a_len = 2.^(len);
    projection = double(0);
    G_inv = zeros(decomp_len);

    % update G array
    for i = 1:decomp_len
        G(i,k) = CH_CH_inner_product(stab_array(i),stab_array(k));
        G(k,i) = CH_CH_inner_product(stab_array(k),stab_array(i));
    end

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
    else % Gram matrix is singular -> find a better stabilizer decomp
        projection = -1;
    end
end
