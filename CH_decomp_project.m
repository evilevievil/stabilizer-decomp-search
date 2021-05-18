%%
%%  projection_array = CH_decomp_project(a,stab_array,len,decomp_len)                                                        
%%  computes projection of a onto stabilizers in the decomp, 
%%  where |a> is arbitrary state, stab_array is array of stabilizers of bit length len
%%
% waiting for CH_CH_inner_product implementation to finish... :(
function projection = CH_decomp_project(a,stab_array,len,decomp_len)
    G = zeros(decomp_len);
    a_len = 2.^(len);
    a_stab_array = zeros(1,decomp_len);

    % compute G array
    for i = 1:decomp_len
        for j=i:decomp_len
            if i==j 
                G(i,j) = 1; % assume unit vector for now
            else
                G(i,j) = CH_CH_inner_product(stab_array(i),stab_array(j));
            end
            G(j,i) = G(i,j);
        end
    end
    G_inv = inv(G);
    
    % compute basis array
    for i = 1:decomp_len
        for j = 1:a_len
            a_stab_array(i) = a_stab_array(i)  + a(j) * CH_basis_inner_product(uint8(j),stab_array(i));
        end
    end

    % compute projection
    for i = 1:decomp_len
        for j = 1:decomp_len
            projection = projection + a_stab_array(i) * a_stab_array(j) * G_inv(i,j);
        end
    end
end