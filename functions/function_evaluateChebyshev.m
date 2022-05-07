function y = function_evaluateChebyshev( mun , x )
    % Evaluates f(x) = (2/pi) / sqrt(1-x^2) * sum_{n=0}^{N-1} mun(n) * T_n(x).
    
    y = zeros(size(x));
    for i = 1 : length(x)
        
        f = 0;
        for n = 1 : length(mun)
            f = f + mun(n)*cos((n-1)*acos(x(i)));
        end
        y(i) = 2/pi / sqrt(1-x(i)^2) * f;
        
    end
    
    return;
end
