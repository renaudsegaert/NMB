clear

f1 = @(x) (x.^10 - 1);
f2 = @(x) (1./(10*x.^2+1));
f3 = @(x) (x.^2);

f = f1;


figure;
hold on
grid on
axis equal 
xlabel('x')
ylabel('y')

fplot(f, [-1, 1], 'r', 'LineWidth', 3);

N = 10;
% conditions = zeros(N, 1);
for n=1:5:61
    disp(n)
    x = linspace(-1, 1, n);
%     x = rootsOfChebyshev(n);

    [c, kappa] = interpolate(x, f);
    
    longerX = linspace(-1, 1, 2*n);

    interpolation = evalCheb2(c, longerX);
    plot(longerX, interpolation, 'g');
    
%     plot(longerX, interpolation - f(longerX)', 'b')


%     [c, kappa] = interpolate(chebyshevRoots, f);
%     conditions(n) = kappa;
%     % Plot the chebyshev-spaced interpolation

end

% plot(1:N, conditions)

title('Plot van chebyshev-nulpunten interpolaties voor f_1.')
legend('Original function', 'Function interpolated')
% title('Error on the n-th Chebyshev-nulpunten interpolation of f_2')
% legend('Error')

