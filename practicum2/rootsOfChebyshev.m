function roots = rootsOfChebyshev(k)
%UNTITLED3 returns a column vector containing the n roots of the n-th
% degree Chebyshev function.

roots = zeros(k,1);
for m=1:k
    roots(m) = cos((2*m+1)*pi/(2*k));
end
end

