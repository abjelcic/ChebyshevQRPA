function [x,y,gamma_smear] = function_readDIRQFAMstrength( path , fileName )
    
    path     = string(path);
    fileName = string(fileName);

    fid = fopen( path + fileName );
    assert( fid ~= -1 , 'Cannot open file: ' + path + fileName );
    assert( extractAfter(fileName,strlength(fileName)-4) == '.out' , 'Strength not in .out format!' );
    
    
    for i = 1 : 8
        Line = fgetl(fid);
    end
    Line        = fgetl(fid);
    Tokens      = split(Line," ");
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
    
