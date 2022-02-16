function [mun,Omegab] = function_readDIRQFAMmun( path , fileName )
    
    fid = fopen( strcat(path,fileName) );
    assert( fid ~= -1 , strcat("Cannot open file: ",path,fileName) );

    
    for i = 1 : 9
        Line = fgetl(fid);
    end
    
    Line   = fgetl(fid);
    Tokens = strsplit(Line," ");
    Omegab = sscanf(Tokens{end-1},"%f");
    
    for i = 1 : 3
        Line = fgetl(fid);
    end
    
    n   = 0;
    mun = [];
    while( true )
        Line = fgetl(fid);
        if( ~ischar(Line) ) break; end
        
        [nmu] = sscanf( Line , "%f %f" );
        nn = nmu(1);
        mu = nmu(2);
        
        assert( nn == n , 'nn!=n' );
        
        mun = [ mun , mu ];
    
        n = n + 1;
    end
    
    
    return;
end
    
