function is_approx_equal = approx_equal(a,b,e)
    result = abs(a-b) <= e + e*(1i);
    is_approx_equal = prod(result);
    %disp(result);disp(is_approx_equal);
end