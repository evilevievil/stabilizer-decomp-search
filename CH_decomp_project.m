%%
%%  projection_array = CH_decomp_project(a,stab_array,len,decomp_len)                                                        
%%  computes projection of a onto stabilizers in the decomp, 
%%  where |a> is arbitrary state, stab_array is array of stabilizers of bit length len
%%
function [projection,G,a_stab_array] = CH_decomp_project(a,stab_array,len,decomp_len)
    G = zeros(decomp_len);
    G_inv = zeros(decomp_len);
    a_len = 2.^(len);
    a_stab_array = zeros(1,decomp_len);
    projection = double(0);

    % compute G array
    for i = 1:decomp_len
        for j=1:decomp_len
            G(i,j) = CH_CH_inner_product(stab_array(i),stab_array(j));
        end
    end

    if abs(det(G)-0) > 0.001
        G_inv = inv(G);
        % compute basis array
        for i = 1:decomp_len
            for j = 1:a_len
                % careful! Matlab 1 indexing :(
                a_stab_array(i) = a_stab_array(i)  + conj(a(j)) * CH_basis_inner_product(j-1,stab_array(i));
            end
        end

        % compute projection
        for i = 1:decomp_len
            for j = 1:decomp_len
                projection = projection + a_stab_array(i) * conj(a_stab_array(j)) * G_inv(i,j);
            end
        end
    else % Gram matrix is singular -> find a better stabilizer decomp
        projection = -1;
    end
end
