function mun = function_applyKernel( mun , kernelName , kernelParameter )
        
    switch( kernelName )
        case 'Lorentz'
            if( exist('kernelParameter','var') )
                lambdaLorentz = kernelParameter;
            else
                lambdaLorentz = 1.0;
            end
            gn = @(n,N) function_gnLorentz(n,N,lambdaLorentz);    
        case 'Jackson'
            gn = @(n,N) function_gnJackson(n,N);
        case 'Dirichlet'
            gn = @(n,N) function_gnDirichlet(n,N);
        case 'Lanczos'
            if( exist('kernelParameter','var') )
                MLanczos = kernelParameter;
            else
                MLanczos = 3.0;
            end
            gn = @(n,N) function_gnLanczos(n,N,MLanczos);
        otherwise
            assert(false,'Unknown kernel');
    end



    N = length(mun);
    for n = 0 : N-1
        mun(n +1) = mun(n +1) * gn(n,N);
    end
    
    return;
end


function res = function_gnDirichlet( n , N )
    res = 1;
    return;
end

function res = function_gnFejer( n , N )
    res = 1-n/N;
    return;
end

function res = function_gnWangZugner( n , N , alpha , beta )
    res = exp( - (alpha * n/N)^beta );
    return;
end

function res = function_gnJackson( n , N )
    res = (N+1-n)*cos(pi*n/(N+1)) + sin(pi*n/(N+1))*cot(pi/(N+1));
    res = res / (N+1);
    return;
end

function res = function_gnLanczos( n , N , M )
    if( n == 0 )
        res = 1;
    else
        res = ( sin(pi*n/N)/(pi*n/N) )^M;
    end
    return;
end

function res = function_gnLorentz( n , N , lambda )
    res = sinh( lambda*(1-n/N) ) / sinh(lambda);
    return;
end


