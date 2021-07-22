%% return true if each a_i, b_i are approx equal within error e + ei
%% used in state vector tests
function is_approx_equal = approx_equal(a,b,e)
    result = abs(a-b) <= e + e*(1i);
    is_approx_equal = prod(result);
end
