function A = function_readmatd( path , fileName )
    
    path     = string(path);
    fileName = string(fileName);

    fid = fopen( path + fileName );
    assert( fid ~= -1 , 'Cannot open file: ' + path + fileName );
    assert( extractAfter(fileName,strlength(fileName)-4) == '.mat' , 'Matrix not in .mat format!' );
    
    
    Line = fgetl(fid);
    Line = fgetl(fid); m = sscanf( Line , "%d" );
    Line = fgetl(fid);
    Line = fgetl(fid); n = sscanf( Line , "%d" );
    
    A = zeros(m,n);
    while( true )
        Line = fgetl(fid);
        if( ~ischar(Line) ) break; end
        
        [ija] = sscanf( Line , "%d %d %f" );
        i = ija(1);
        j = ija(2);
        a = ija(3);
        
        A(i,j) = a;
    end
    
    
    return;
end
