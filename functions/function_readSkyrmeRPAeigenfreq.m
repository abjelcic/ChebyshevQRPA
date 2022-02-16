function [Omegas,BIS,BIV] = function_readSkyrmeRPAeigenfreq( path , fileName )
    
    fid = fopen( strcat(path,fileName) );
    assert( fid ~= -1 , strcat("Cannot open file: ",path,fileName) );
    
    
    while( true )
        Line = fgetl(fid);
        if( ~ischar(Line) )
            assert( false , 'Error in function_readSkyrmeRPAeigenfreq' );
        end
        if( strcmp( Line , '          E          B(IS)       %M0(IS)       B(EM)        B(IV)       %M0(IV)' ) == 1 )
            Line = fgetl(fid);
            Line = fgetl(fid);
            break;
        end
    end
    
    
    Omegas = [];
    BIS    = [];
    BIV    = [];
    while( true )
        Line = fgetl(fid);
        if( ~ischar(Line) ) break; end
        if( length(Line(find(~isspace(Line))))==0 ) break; end
        
        [line] = sscanf( Line , "%d %f %f %f %f %f %f" );
        omega = line(2);
        bis   = line(3);
        biv   = line(6);
                
        Omegas = [ Omegas , omega ];
        BIS    = [ BIS , bis ];
        BIV    = [ BIV , biv ];
    end
    
    return;
end
    
