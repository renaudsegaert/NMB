function v = evalCheb(a,x)
    [~,n] = size(a);
    [~,nx] = size(x);
    v = zeros(nx,1);
    for j = 1:nx
        T = zeros(n,1);
        som = 0;
        syms y
        for i = 1:n
            if i == 1
                T(i) = 1 ; 
                som =  som + a(i)*T(i);
            elseif i ==2
                T(i)  = x(j); 
                som =  som + a(i)*T(i);
            else
                
                T(i) = 2*x(j)*T(i-1) - T(i-2);
                som = som + a(i)*T(i);
            end

        end
        
        v(j) = som;
    end
    
    
    
    
    
end