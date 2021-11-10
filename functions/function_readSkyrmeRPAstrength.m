function [x,y] = function_readSkyrmeRPAstrength( path , fileName )
    
    path     = string(path);
    fileName = string(fileName);

    fid = fopen( path + fileName );
    assert( fid ~= -1 , 'Cannot open file: ' + path + fileName );
    assert( extractAfter(fileName,strlength(fileName)-4) == '.dat' , 'Strength not in .dat format!' );
    
    x = [];
    y = [];
    while( true )
        Line = fgetl(fid);
        if( ~ischar(Line) ) break; end
        
        [wdBdw] = sscanf( Line , "%f %f" );
        w    = wdBdw(1);
        dBdw = wdBdw(2);
                
        x = [ x , w    ];
        y = [ y , dBdw ];
    end
    
    return;
end
    
