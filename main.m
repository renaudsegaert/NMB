%bronnen:
% gram schmidt https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/768d6fe735722547c3f5e74824532c53_MIT18_06S10_gramschmidtmat.pdf
% gram schmidt https://www.math.ucla.edu/~yanovsky/Teaching/Math151B/handouts/GramSchmidt.pdf
% https://arnold.hosted.uark.edu/NLA/Pages/CGSMGS.pdf
% https://math.stackexchange.com/questions/4425346/how-to-prove-that-the-sum-of-two-random-high-dimensional-vectors-is-close-to-s

%% Opgave 1

%OPDRACHT 3
B = randn(200,100);
[Q,~] = qr(B,0);
D = zeros(100);
for i=1:100
   D(i,i) = 2^(1-i);
end
V = zeros(100,100);
for i = 1:100
    V(i,i) = 1;
    V(1,i)=1;    
end
A = Q*D*V;


%OPDRACHT 4
[Qkl , Rkl] = klGS(A);
[Qmod , Rmod] = modGS(A);
a = 70; 
Diffkl = zeros(a,1);
for i = 1:a
    Diffkl(i,1) = norm(Q(:,i)-Qkl(:,i))  ;
end
Diffmod = zeros(a,1);
for i = 1:a
    Diffmod(i,1) = norm(Q(:,i)-Qmod(:,i))  ;
end
xas = (1:1:a);
figure
p = semilogy(xas,Diffkl,xas,Diffmod);
title('Error between Q and calculated Q from algorithms');
xlabel('$q_i$','Interpreter','Latex') 
ylabel('$ |\tilde{q_i}-q_i|_2$','Interpreter','Latex') 







%OPDRACHT 5
Difftheoriemod = zeros(a,1);
for i = 1:a
    b = eps*2^(i-1)*Q(:,1);
    Difftheoriemod(i,1) =  norm((Q(:,i)-b)-Q(:,i));
end
Difftheoriekl = zeros(a,1);

for i = 1:a
    summatie = 0;
    for j = 2:i-1
        summatie = summatie + Q(:,j)*eps*2^(j-1); 
    end
    
    Difftheoriekl(i,1) = norm( -1*(2^(i-1)*eps*Q(:,1))+ 2^(i-1)*summatie);     
end
xas = (1:1:a);
figure
p = semilogy(xas,Diffkl,xas,Diffmod,xas,Difftheoriekl,xas,Difftheoriemod);
title('Error Q and ~Q in theory and in practice');
xlabel('$q_i$','Interpreter','Latex') 
ylabel('$ |\tilde{q_i}-q_i|_2$','Interpreter','Latex')




%OPDRACHT 7
Difftheoriemodecht = zeros(a,1);
for i = 1:a
    b = norm(((Q(:,i)/2^(i-1))-eps*Q(:,1)));
    c = ((Q(:,i)/2^(i-1))-eps*Q(:,1));
    Z = c/b;
    Difftheoriemodecht(i,1) = norm(Q(:,i)-Z);    
end
xas = (1:1:a);
figure
p = semilogy(xas,Diffmod,xas,Difftheoriemodecht);
title('Error Q and ~Q in theory and in practice for big n');
xlabel('$q_i$','Interpreter','Latex') 
ylabel('$ |\tilde{q_i}-q_i|_2$','Interpreter','Latex')



%OPDRACHT 8
Difforthogonaliteitkl = zeros(a,1);
for i = 1:a
    Difforthogonaliteitkl(i,1) = norm(eye(i)-(Qkl(:,1:i)'*Qkl(:,1:i)));
end
Difforthogonaliteitmod = zeros(a,1);
for i = 1:a
    Difforthogonaliteitmod(i,1) = norm(eye(i)-(Qmod(:,1:i)'*Qmod(:,1:i)));
end
xas = (1:1:a);
figure
p = semilogy(xas,Diffkl,xas,Diffmod,xas,Difforthogonaliteitkl,xas,Difforthogonaliteitmod);
title('Loss of orthogonality from classical gram-schmidt and modified gram-schmidt');
xlabel('$q_i$','Interpreter','Latex') 
ylabel('$ |I-\tilde{Q_i}^{*}\tilde{Q_i}|_2$','Interpreter','Latex')



%OPDRACHT 9
B = randn(200,50);
C = randn(50,50);
[V,~] = qr(C);
[Q,~] = qr(B,0);
D = diag(2.^linspace(0,1,50));
A = Q*D*V;
disp(cond(A));



%OPDRACHT 10
Diffkkl = zeros(50,1);
Diffkmod = zeros(50,1);
Diffkher = zeros(50,1);
B = randn(200,50);
[Q,R] = qr(B,0);
C = randn(50,50);
[V,~] = qr(C);
for i = 1:50
    D = diag(2.^linspace(0,i,50));
    A = Q*D*V;    
    [Qkl,Rkl] = klGS(A);
    [Qmod,Rmod] = modGS(A);
    [Qher,Rher] = herGS(A);
    Diffkkl(i,1) = norm(eye(50)-(Qkl'*Qkl));
    Diffkmod(i,1) = norm(eye(50)-(Qmod'*Qmod)); 
    Diffkher(i,1) = norm(eye(50)-(Qher'*Qher));
end
xas = (1:1:50);
figure
semilogy(xas,Diffkkl,xas,Diffkmod,xas,Diffkher)
title('Loss of orthogonality depending on k','Interpreter','Latex')
xlabel('k','Interpreter','Latex') 
ylabel('$|I-Q^{*}Q|$','Interpreter','Latex')




%OPDRACHT 12
%
%



%% Opgave 2
% TODO: Deal with divergence of the maximum eigen value.

clear
close all
hold on

nb_eigenvalues_to_calculate = 6;
k = 1000;
A = sprand(k,k, 1/k);
maxA = max(eigs(A, nb_eigenvalues_to_calculate));
v = rand(k,1);

% plot(real(maxA), imag(maxA), 'gO');

for m=1:100
    disp(m);
    [V, H] = Arnoldi(A,v,m);
    sh = size(H);
    
    eigenvalues_H = eigs(H(1:m,:));
    for j=1:length(eigenvalues_H)
        current_eigenvalue = eigenvalues_H(j);
        plot(m, real(current_eigenvalue), 'g*');
    end
end

eigenvalues_A = eigs(A);
for j=1:length(eigenvalues_A)
    current_eigenvalue = eigenvalues_A(j);
    plot(m+1, real(current_eigenvalue), 'rO');
end

eigs(A)
eigs(H(1:sh(1)-1,:))

hold off
title("Ritzwaarden per Iteratiestap");
xlabel('Iteratiestap');
ylabel('Real');





