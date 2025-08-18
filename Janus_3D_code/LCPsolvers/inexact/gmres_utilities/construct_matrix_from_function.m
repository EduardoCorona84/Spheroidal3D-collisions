function operator = construct_matrix_from_function(f, input_size)
    %Constructs a matrix representation of the linear operator f
    identity = eye(input_size);
    firstColumn = f(identity(:,1));
    output_size = size(firstColumn, 1);
    operator = zeros(output_size, input_size);
    operator(:,1) = firstColumn;
    for i = 2:input_size
        operator(:,i) = f(identity(:,i)); 
    end

end