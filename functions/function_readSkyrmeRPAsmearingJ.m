function [gamma,J] = function_readSkyrmeRPAsmearingJ( path , fileName )
    
    path     = string(path);
    fileName = string(fileName);

    fid = fopen( path + fileName );
    assert( fid ~= -1 , 'Cannot open file: ' + path + fileName );
    assert( extractAfter(fileName,strlength(fileName)-3) == '.in' , 'Smearing not in .in format!' );
    
   
    for i = 1 : 6
        Line = fgetl(fid); 
    end
    
    Line = fgetl(fid);
    [line] = sscanf( Line , "%d %d" );
    J = line(1);
    
    for i = 1 : 2
        Line = fgetl(fid); 
    end  
    
    Line = fgetl(fid);
    [line] = sscanf( Line , "%f %f %f %f" );
    
    Gamma = line(4);
    gamma = 0.5 * Gamma;
    
    return;
end
    
