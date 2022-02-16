function [x,y] = function_readSkyrmeRPAstrength( path , fileName )
    
    fid = fopen( strcat(path,fileName) );
    assert( fid ~= -1 , strcat("Cannot open file: ",path,fileName) );
    
    
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
    
