function V = function_readvecd( path , fileName )
    
    fid = fopen( strcat(path,fileName) );
    assert( fid ~= -1 , strcat("Cannot open file: ",path,fileName) );
    
    
    Line = fgetl(fid);
    Line = fgetl(fid); n = sscanf( Line , "%d" );
    
    V = zeros(n,1);
    while( true )
        Line = fgetl(fid);
        if( ~ischar(Line) ) break; end
        
        [iv] = sscanf( Line , "%d %f" );
        i = iv(1);
        v = iv(2);
        
        V(i) = v;
    end
    
    return;
end