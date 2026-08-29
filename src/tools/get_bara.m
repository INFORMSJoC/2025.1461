function bara = get_bara(A,Aidx)

if ~iscell(A)
    Atemp{1} = A;
    A = Atemp;
end

if nargin > 1

    % Initialize empty arrays with an estimated size
    [ii, jj] = size(A);
    nnz_est = 0;
    for i = 1:ii
        for j = 1:jj
            if Aidx(i,j)
                nnz_est = nnz_est + nnz(tril(A{i,j}));
            end
        end
    end

else


    % Initialize empty arrays with an estimated size
    [ii, jj] = size(A);
    nnz_est = 0;
    for i = 1:ii
        for j = 1:jj
            if ~isempty(A{i,j})
                nnz_est = nnz_est + nnz(tril(A{i,j}));
            end
        end
    end

end

% tic
% [row_idx, col_idx] = find(~cellfun(@isempty, As));
% for k = 1:length(row_idx)
%     A = As{row_idx(k), col_idx(k)};
%     nnz_est = nnz_est + nnz(tril(A));
% end
%
% toc

% Preallocate arrays
bara.subi = zeros(1, nnz_est);
bara.subj = zeros(1, nnz_est);
bara.subk = zeros(1, nnz_est);
bara.subl = zeros(1, nnz_est);
bara.val  = zeros(1, nnz_est);

if nargin > 1
    % Fill arrays efficiently
    idx = 1;
    for i = 1:ii
        for j = 1:jj
            if Aidx(i,j)
                [rowA, colA, valA] = find(tril(A{i,j}));  % Extract lower-triangle


                m = length(rowA);
                if m > 0
                    idx_end = idx + m - 1;

                    bara.subi(idx:idx_end) = i;  % Constraint index
                    bara.subj(idx:idx_end) = j;  % Constraint index
                    bara.subk(idx:idx_end) = rowA;  % Row index
                    bara.subl(idx:idx_end) = colA;  % Column index
                    bara.val(idx:idx_end) = valA;  % Value

                    idx = idx_end + 1;
                end
            end
        end
    end

else

    % Fill arrays normally
    idx = 1;
    for i = 1:ii
        for j = 1:jj
            if ~isempty(A{i,j})
                [rowA, colA, valA] = find(tril(A{i,j}));  % Extract lower-triangle


                m = length(rowA);
                if m > 0
                    idx_end = idx + m - 1;

                    bara.subi(idx:idx_end) = i;  % Constraint index
                    bara.subj(idx:idx_end) = j;  % Constraint index
                    bara.subk(idx:idx_end) = rowA;  % Row index
                    bara.subl(idx:idx_end) = colA;  % Column index
                    bara.val(idx:idx_end) = valA;  % Value

                    idx = idx_end + 1;
                end
            end
        end
    end

end

% Trim unused preallocated space
bara.subi = bara.subi(1:idx-1);
bara.subj = bara.subj(1:idx-1);
bara.subk = bara.subk(1:idx-1);
bara.subl = bara.subl(1:idx-1);
bara.val  = bara.val(1:idx-1);

end