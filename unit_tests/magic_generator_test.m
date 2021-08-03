



gen_array = zeros(6,1);
gen_array = cast(gen_array,const.typecast_str);
gen_array(5) = cast(bin2dec('0000000000010001'),const.typecast_str);
gen_array(3) = cast(bin2dec('0000000000000101'),const.typecast_str);
leading_bits = cast(bin2dec('0000000000010100'),const.typecast_str);

x_bits = get_commuter(gen_array,leading_bits,6);

dec2bin(x_bits)


