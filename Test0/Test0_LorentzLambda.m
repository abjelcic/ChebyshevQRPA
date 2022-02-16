% MATLAB R2018a
clear;
clc;
close all;
addpath('../functions');


figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,1,1.0,1.0]);


NPlotPoints = [ 1000 , 500 , 250 , 100 , 50 ];

lambdas = [ 1.0 , 2.0 , 3.0 , 4.0 , 5.0 ];
for i = 1 : length(lambdas)
    
    lambda = lambdas(i);
    
    N = 64;
    mun = cos( pi/2 * [0:N-1] ) ./ ( 1 + ([0:N-1]==0) );
    mun = function_applyKernel( mun , 'Lorentz' , lambda );
    
    
    gamma = lambda/N;
    xPlotMesh1 = linspace(-0.5,+0.5,NPlotPoints(i));
    yPlotMesh1 = (gamma/pi) ./ ( gamma^2 + xPlotMesh1.^2 );
    
    [xPlotMesh2,yPlotMesh2] = function_fftEvaluateChebSeries( mun , 1.0 , 100000 );
    
 
    subplot(1,length(lambdas),i);
    plot( xPlotMesh1 , yPlotMesh1 , 'ro' , 'LineWidth' , 2.0 , 'MarkerSize' , 3 ); hold on;
    plot( xPlotMesh2 , yPlotMesh2 , 'b-' , 'LineWidth' , 1.0 , 'MarkerSize' , 3 ); hold on;
    
    
    xlim([-0.5,+0.5]);
    xlabel({"$x$",strcat("$\lambda = ",num2str(lambda)," $")},'Interpreter','latex');
    xticks([ -0.50 , 0.00 , +0.50 ]);

    ylim([0,20]);
    if( i == 1 ) ylabel('$y$','Interpreter','latex'); end
    
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'XMinorTick','on','YMinorTick','on');
        
end
