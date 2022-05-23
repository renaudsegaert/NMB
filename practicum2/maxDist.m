function d = maxDist(f1,f2)
%Bepaalt de maximale afstand tussen twee functies in het interval [-1, 1]

x = linspace(-1, 1, 200)';
differences = abs(f1(x) - f2(x));
d = max(differences);
end

