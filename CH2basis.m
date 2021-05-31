function basis_vector = CH2basis(stab)
    vec_len = 2.^(stab.len);
    basis_vector = zeros(vec_len,1);
    for i = 1:vec_len
        %% matlab 1 indexing!! :(
        %% reverse bit string due right most convention
        diff_len = stab.len - strlength(dec2bin(i-1));
        i_strrev = reverse(dec2bin(i-1));
        for j = 1: diff_len
            i_strrev = append(i_strrev,'0');
        end
        basis_vector(i) = CH_basis_inner_product(bin2dec(i_strrev),stab);
    end
end