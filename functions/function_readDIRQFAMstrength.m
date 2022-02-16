function [x,y,gamma_smear] = function_readDIRQFAMstrength( path , fileName )
    
    fid = fopen( strcat(path,fileName) );
    assert( fid ~= -1 , strcat("Cannot open file: ",path,fileName) );
    
    
    for i = 1 : 8
        Line = fgetl(fid);
    end
    Line        = fgetl(fid);
    Tokens      = strsplit(Line," ");
    gamma_smear = sscanf(Tokens{end-1},"%f");
    for i = 1 : 4
        Line = fgetl(fid);
    end    
    
    
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
    
