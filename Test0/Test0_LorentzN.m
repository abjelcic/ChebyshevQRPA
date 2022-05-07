clear;
clc;
close all;
addpath('../functions');


figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,1,1.0,1.0]);


NPlotPoints = [ 50 , 100 , 250 , 500 , 1000 ];

Nsizes = [ 16 , 32 , 64 , 128 , 256 ];
for i = 1 : length(Nsizes)
    
    N = Nsizes(i);
    
    lambda = 1.5;
    mun = cos( pi/2 * [0:N-1] ) ./ ( 1 + ([0:N-1]==0) );
    mun = function_applyKernel( mun , 'Lorentz' , lambda );
    
    
    gamma = lambda/N;
    xPlotMesh = linspace(-0.5,+0.5,NPlotPoints(i));
    
    yPlotMesh1 = (gamma/pi) ./ ( gamma^2 + xPlotMesh.^2 );
    yPlotMesh2 = function_evaluateChebyshev( mun , xPlotMesh );
    
 
    subplot(1,length(Nsizes),i);
    plot( xPlotMesh , yPlotMesh1 , 'ro' , 'LineWidth' , 2.0 , 'MarkerSize' , 3 ); hold on;
    plot( xPlotMesh , yPlotMesh2 , 'b-' , 'LineWidth' , 1.0 , 'MarkerSize' , 3 ); hold on;
    
    
    xlim([-0.5,+0.5]);
    xlabel({"$x$",strcat("$N = ",num2str(N),"$")},'Interpreter','latex');
    xticks([ -0.50 , 0.00 , +0.50 ]);

    ylim([0,55]);
    if( i == 1 ) ylabel('$y$','Interpreter','latex'); end

    set(gca,'TickLabelInterpreter','latex');
    set(gca,'XMinorTick','on','YMinorTick','on');

end






function y = function_evaluateChebyshev( mun , x )
    % Evaluates f(x) = (2/pi) / sqrt(1-x^2) * sum_{n=0}^{N-1} mun(n) * T_n(x).
    
    y = zeros(size(x));
    for i = 1 : length(x)
        
        f = 0;
        for n = 1 : length(mun)
            f = f + mun(n)*cos((n-1)*acos(x(i)));
        end
        y(i) = 2/pi / sqrt(1-x(i)^2) * f;
        
    end
    
    return;
end
