function tensor_prod = tensor_exp(A,n)
    tensor_prod = 1;
    for i = 1:n
        tensor_prod = kron(tensor_prod,A);
    end
end