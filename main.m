%bronnen:
% gram schmidt https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/768d6fe735722547c3f5e74824532c53_MIT18_06S10_gramschmidtmat.pdf
% gram schmidt https://www.math.ucla.edu/~yanovsky/Teaching/Math151B/handouts/GramSchmidt.pdf
% https://arnold.hosted.uark.edu/NLA/Pages/CGSMGS.pdf
%https://math.stackexchange.com/questions/4425346/how-to-prove-that-the-sum-of-two-random-high-dimensional-vectors-is-close-to-s


%OPDRACHT 3
B = randn(200,100);

[Q,R] = qr(B,0);

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

% Antwoord: stel A = Q*D*V --> a1 = q1,  ai = (qi*2^(1-i))+q1

%OPDRACHT 4

[Qkl , Rkl] = klGS(A);
[Qmod , Rmod] = modGS(A);
[Qher , Rher] = herGS(A);
a = 100; 
Diffkl = zeros(a,1);
for i = 1:a
    Diffkl(i,1) = norm(Q(:,i)-Qkl(:,i))  ;
end
Diffmod = zeros(a,1);
for i = 1:a
    Diffmod(i,1) = norm(Q(:,i)-Qmod(:,i))  ;
end
Diffher = zeros(a,1);
for i = 1:a
    Diffher(i,1) = norm(Q(:,i)-Qher(:,i))  ;
end


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
Difftheoriemodecht = zeros(a,1);

for i = 1:a
    b = norm(((Q(:,i)/2^(i-1))-eps*Q(:,1)));
    c = ((Q(:,i)/2^(i-1))-eps*Q(:,1));
    Z = c/b;
    Difftheoriemodecht(i,1) = norm(Q(:,i)-Z);    
end

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


plot(xas,Diffkl,xas,Diffmod,xas,Diffher)
title('Error between Q and ~Q')
ylim([-0.2 1.8])
% Antwoord: DE stijging is exponentieel, fout tussen Q en ~Q stijgt eerst
% zachtjes maar met elke stap stijgt hij sneller en sneller tot het het
% punt van afvlakking heeft bereikt.

%OPDRACHT 5
% foutenanalyse klassieke Gramm schmidt:

% ~q1 ? q1

% ã2 ? a2 - <a2,~q1>*~q1
% ã2 ? q1 + q2/2 -(1+?1)*q1
% ã2 = q2/2 -?1*q1
% ~q2 = q2 - 2*?1*q1

%ã3 = a3 - <a3,~q1>*~q1 - <a3,~q2>*~q2
%ã3 = q1 + q3/4 - (1 + ?2)*q1 - (-2*?1)*q2
%ã3 = q3/4 - ?2*q1 + 2*?1*q2
%~q3 = q3 - 4*?2*q1 + 8*?1*q2

%~q4 = q4 - 8*?1*q1 + 16*?2*q2+24*?3*q3
% voor n > 2 
%~qn = qn +       2^(n-1)*((-?1*q1)+summation(from i=2 to i =n (2^(n-1)*?i*q1))

% foutenanalyse aangepaste Gramm schmidt:

%~qn = qn - (2^(n-1)*?n*q1)


%OPDRACHT 6
%

%OPDRACHT 7
% Voor grote n is (qn/2^(n-1))-?n*q1 niet meer ongeveer gelijk aan
% (qn/2^(n-1)), Voor de norm van ãn mogen we dus niet zomaar 1/2^(n-1) pakken

%~qn = ((qn/2^(n-1))-?n*q1)/norm(((qn/2^(n-1))-?n*q1))
% Deze schatting komt in de buurt van de experimentele waarden, beide
% stijgen exponentieel rond het zelfde punt en vlakken af naar het einde
% toe

%OPDRACHT 8
%Wanneer Qj perfect orthogonaal is dan zou Qj'*Qj gelijk moeten zijn aan de
%eenheidsmatrix, om het verschil tussen orthogonaliteit te meten bekijken
%we dus het verschil tss Qj'*Qj en de eenheidsmatrix

%OPDRACHT 9
B = randn(200,50);

[Q,R] = qr(B,0);

D = diag(2.^linspace(0,3,50));
V = eye(50);
A = Q*D*V;
[Qkl,Rkl] = klGS(A);
[Qmod,Rmod] = modGS(A);
[Qher,Rher] = herGS(A);
%condition of a matrix is the highest singular value divided by the lowest
%singular value, the highest value is 2^k and the lowest is always 1 so the
%condition is equal to 2^k

%OPDRACHT 10
Diffkkl = zeros(50,1);
Diffkmod = zeros(50,1);
Diffkher = zeros(50,1);

for i = 1:50
    [Q,R] = qr(B,0);

    D = diag(2.^linspace(0,i,50));
    A = Q*D*V;
    disp(cond(A))
    [Qkl,Rkl] = klGS(A);
    [Qmod,Rmod] = modGS(A);
    [Qher,Rher] = herGS(A);
    Diffkkl(i,1) = norm(eye(50)-(Qkl'*Qkl));
    Diffkmod(i,1) = norm(eye(50)-(Qmod'*Qmod)); 
    Diffkher(i,1) = norm(eye(50)-(Qher'*Qher));
end


xas = (1:1:50);
figure
plot(xas,Diffkkl,xas,Diffkmod,xas,Diffkher)
title('norm(I-Q(transposed)*Q')



