function V = function_readvecd( path , fileName )
    
    path     = string(path);
    fileName = string(fileName);

    fid = fopen( path + fileName );
    assert( fid ~= -1 , 'Cannot open file: ' + path + fileName );
    assert( extractAfter(fileName,strlength(fileName)-4) == '.vec' , 'Vector not in .vec format!' );
    
    
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