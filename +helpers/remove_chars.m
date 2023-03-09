function C = remove_chars(C, n)
% Removes first n characters from every string element
% belonging to the 1 x N cell C
% Not vectorized
for idx = 1:size(C,2)
    C{idx} = C{idx}((n+1):end);
end
end