function matrices = construct_test_matrices(problem_size)
    %This function constructs a set of 4 test matrices for the LCP solver of a specific size. They are (based off xtrace spectrum paper):
    % 1. Linear condition 50: lambda = 50 - (i - 1)/(N - 1) for i = 1:problem_size
    % 2. Linear condition size: lambda = problem_size - (i - 1)/(N - 1) for i = 1:problem_size
    % 3. Flat: lambda = 3 - 2(i - 1)/(N - 1) for i = 1:problem_size
    % 4. exp:  lambda = exp(0.7 ^(i - 1)) for i = 1:problem_size

    %linear condition 50
    matrices.linear50.eigenvalues = linspace(50, 1, problem_size);
    matrices.linear50.eigenvectors = orth(randn(problem_size));
    matrices.linear50.matrix = matrices.linear50.eigenvectors * diag(matrices.linear50.eigenvalues) * matrices.linear50.eigenvectors';
    matrices.linear50.inverse = matrices.linear50.eigenvectors * diag(1./matrices.linear50.eigenvalues) * matrices.linear50.eigenvectors';
    matrices.linear50.matrixNorm = max(matrices.linear50.eigenvalues);
    matrices.linear50.inverseNorm = 1/min(matrices.linear50.eigenvalues);

    %linear condition size
    matrices.linearSize.eigenvalues = linspace(problem_size, 1, problem_size);
    matrices.linearSize.eigenvectors = orth(randn(problem_size));
    matrices.linearSize.matrix = matrices.linearSize.eigenvectors * diag(matrices.linearSize.eigenvalues) * matrices.linearSize.eigenvectors';
    matrices.linearSize.inverse = matrices.linearSize.eigenvectors * diag(1./matrices.linearSize.eigenvalues) * matrices.linearSize.eigenvectors';
    matrices.linearSize.matrixNorm = max(matrices.linearSize.eigenvalues);
    matrices.linearSize.inverseNorm = 1/min(matrices.linearSize.eigenvalues);
    
    %flat
    matrices.flat.eigenvalues = zeros(problem_size, 1);
    for i = 1:problem_size
        matrices.flat.eigenvalues(i) = 3 - 2*(i - 1)/(problem_size - 1);
    end
    matrices.flat.eigenvectors = orth(randn(problem_size));
    matrices.flat.matrix = matrices.flat.eigenvectors * diag(matrices.flat.eigenvalues) * matrices.flat.eigenvectors';
    matrices.flat.inverse = matrices.flat.eigenvectors * diag(1./matrices.flat.eigenvalues) * matrices.flat.eigenvectors';
    matrices.flat.matrixNorm = max(matrices.flat.eigenvalues);
    matrices.flat.inverseNorm = 1/min(matrices.flat.eigenvalues);

    %exponential
    matrices.exp.eigenvalues = zeros(problem_size, 1);
    for i = 1:problem_size
        matrices.exp.eigenvalues(i) = exp(0.7^(i - 1));
    end
    matrices.exp.eigenvectors = orth(randn(problem_size));
    matrices.exp.matrix = matrices.exp.eigenvectors * diag(matrices.exp.eigenvalues) * matrices.exp.eigenvectors';
    matrices.exp.inverse = matrices.exp.eigenvectors * diag(1./matrices.exp.eigenvalues) * matrices.exp.eigenvectors';
    matrices.exp.matrixNorm = max(matrices.exp.eigenvalues);
    matrices.exp.inverseNorm = 1/min(matrices.exp.eigenvalues);


end