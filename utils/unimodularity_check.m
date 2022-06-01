function [is_unimodular] = unimodularity_check(A)

    %A = dlmread('results/numeric_results/constraints_matrix.txt');
    %A = A(2:end, :);
    disp(size(A));
    [rows, cols] = size(A);
    % disp(rank(A));
    is_unimodular = true;

    sprank_A = sprank(A);
    if sprank_A == cols
        disp('sprank = cols');
        rank_A = sprank_A;
    else
        rank_A = cols;
    end
    
    for start = 1: rows + 1 - rank_A
        
        if mod(start, 10000) == 0
            disp('ind');
            disp(start);
        end
        
        if sprank(A(start : start + rank_A - 1, :)) < cols
            continue;
        end
        d = det(A(start : start + rank_A - 1, :));
        if d ~= 0 && d ~= 1 && d ~= -1
            disp('det');
            disp(d);
            disp('start, end');
            disp([start, start+rank_A-1]);
            is_unimodular = false;
        end
    end
end