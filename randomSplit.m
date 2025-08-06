function [set1, set2] = randomSplit(n)
    % Generate IDs from 1 to n
    ids = randperm(n);  % Random permutation of integers from 1 to n
    
    % Determine split point
    splitPoint = ceil(n / 2);
    
    % Split the randomized IDs
    set1 = ids(1:splitPoint);
    set2 = ids(splitPoint+1:end);
end
