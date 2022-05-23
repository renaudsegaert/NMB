clear
tic

f1 = @(x) (x.^10 - 1);
f2 = @(x) (1./(10*x.^2+1));
f3 = @(x) (x.^2);

f = f2;


figure;
hold on
grid on
axis equal 
xlabel('Graad n')
ylabel('Maximale fout')
% fplot(f, [-1, 1], 'r', 'LineWidth', 3);

N = 101;
skip = 10;
% conditions1 = zeros(N, 1);
% conditions2 = zeros(N, 1);
maxError = zeros(N, 1);
for n=1:skip:N
    disp(n)
    x = linspace(-1, 1, n);
%     x = rootsOfChebyshev(n);

    [c, kappa] = interpolate(x, f);
    func = @(t) (evalCheb2(c, t));
    maxError(n) = maxDist(f, func);
end

plot(1:skip:N, maxError(1:skip:N) ./ f2(1:skip:N) , 'ro')

% semilogy(1:N, conditions1, 'g*')
% semilogy(1:N, conditions2, 'ro')
% title('Condities van M-matrices van Chebyshev-nulpunten voor interpolatie.')
% % legend('Equidistant', 'Chebyshev')
% xlabel('Graad n')
% ylabel('Conditiegetal')

% title('Plot van chebyshev-nulpunten interpolaties voor f_1.')
% legend('Original function', 'Function interpolated')
% title('Error on the n-th equidistant interpolation of f_1')
% legend('Error')

