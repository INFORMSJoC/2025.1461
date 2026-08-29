function K = lower_tri_index(I, J, n)
    % Computes the position k in the vectorized lower-triangular part
    % of an n x n matrix (column-major order) for entry (i,j), where i >= j

    if any(I < J)
        error('Only lower triangular entries are allowed (i >= j).');
    end

    K = (J-1).*(2*n-J+2)/2 + (I - J + 1);
end
