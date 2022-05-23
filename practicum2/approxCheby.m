function a = approxCheby(f,n)
v = zeros(1,n+1);
b = zeros(1,n+1);
for i = 1:n+1
    b(1,i) = cos((i-1)*pi/n);
    v(1,i) = f( cos((i-1)*pi/n));
end   
w = fliplr(v(2:end-1));
v_even = [v w];

V = fft(v_even)/n;
V = real(V(1:n+1));
V = [V(1)/2 V(2:end)];
a = V;



end