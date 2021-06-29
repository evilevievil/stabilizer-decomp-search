%%
%%  projection_array = CH_decomp_project(a,stab_array,len,decomp_len)                                                        
%%  computes projection of a onto stabilizers in the decomp, 
%%  where |a> is arbitrary state, stab_array is array of stabilizers of bit length len
%%
function projection = CH_decomp_project(a,stab_array,G,G_inv,a_stab_array,projection,k,len,decomp_len)
    %G = zeros(decomp_len);
    %a_stab_array = zeros(1,decomp_len);
    a_len = 2.^(len);
    projection = double(0);

    % update projection
    for j = 1:decomp_len
        projection = projection - a_stab_array(k) * conj(a_stab_array(j)) * G_inv(k,j);
        projection = projection - a_stab_array(j) * conj(a_stab_array(k)) * G_inv(j,k);
    end

    % update G array
    for i = 1:decomp_len
        % todo: just use conjugate for lower triangular
        G(i,k) = CH_CH_inner_product(stab_array(i),stab_array(k));
        G(k,i) = CH_CH_inner_product(stab_array(k),stab_array(i));
    end
    % todo: can I optimize this?
    if abs(det(G)-0) > 0.001
        G_inv = inv(G);
        % update basis array
        for j = 1:a_len
            % careful! Matlab 1 indexing :(
            a_stab_array(k) = a_stab_array(k)  + conj(a(j)) * CH_basis_inner_product(j-1,stab_array(k));
            %disp(a_stab_array);
        end
        
        %disp(a_stab_array);
        %disp(G);

        % update projection
        for j = 1:decomp_len
            projection = projection + a_stab_array(k) * conj(a_stab_array(j)) * G_inv(k,j);
            projection = projection + a_stab_array(j) * conj(a_stab_array(k)) * G_inv(j,k);
        end
    else
        % Gram matrix is singular -> find a better stabilizer decomp
        projection = -1;
    end
end
