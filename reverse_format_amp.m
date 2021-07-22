%% sorts user input amplitude array into reversed bit string order
%% caveat: CH_state indexes bits from right to left, reverse_format_amp 
function reverse_formatted_a = reverse_format_amp(a,len)
    vec_len = 2.^len;
    reverse_formatted_a = zeros(vec_len,1);
    for i = 1:vec_len
        diff_len = len - strlength(dec2bin(i-1));
        i_strrev = reverse(dec2bin(i-1));
        for j = 1: diff_len
            i_strrev = append(i_strrev,'0');
        end
        reverse_formatted_a(i,1) = a(bin2dec(i_strrev)+1,1);
    end
end