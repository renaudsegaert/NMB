function prod = evalCheb2(a, x)
    na = length(a);
    nx = length(x);
    T = zeros(nx, na);

   for j=1:na
       T(:, j) = chebyshevT(j-1, x);
   end
    prod = T * a;
end

