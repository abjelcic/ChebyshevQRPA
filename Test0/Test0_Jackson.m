clear;
clc;
close all;
addpath('../functions');


figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,1,1.0,1.0]);


    
N = 64;
mun = cos( pi/2 * [0:N-1] ) ./ ( 1 + ([0:N-1]==0) );

munJackson   = function_applyKernel( mun , 'Jackson'   );
munDirichlet = function_applyKernel( mun , 'Dirichlet' );


sigma = pi/N;
xPlotMesh  = linspace(-1,+1,300);

yPlotMesh1 = 1/sqrt(2*pi*sigma^2) * exp( -0.5 * ( xPlotMesh./sigma ).^2 );
yPlotMesh2 = function_evaluateChebyshev( munJackson   , xPlotMesh );
yPlotMesh3 = function_evaluateChebyshev( munDirichlet , xPlotMesh );


plot( xPlotMesh , yPlotMesh1 , 'ro'  , 'LineWidth' , 2.0 , 'MarkerSize' , 3 ); hold on;
plot( xPlotMesh , yPlotMesh2 , 'b-'  , 'LineWidth' , 1.0 , 'MarkerSize' , 3 ); hold on;
plot( xPlotMesh , yPlotMesh3 , 'k--' , 'LineWidth' , 1.0 , 'MarkerSize' , 3 ); hold on;


legend1 = "Gaussian $\sigma = \pi/N$";
legend2 = strcat( "KPM ($N=" , num2str(N) , "$)" );
legend3 = "Without KPM";
legend({legend1,legend2,legend3},'Interpreter','latex');
    

xlim([-1,+1]);
xlabel('$x$','Interpreter','latex');
xticks([ -1.00 , -0.75 , -0.50 , -0.25 , 0.00 , +0.25 , +0.50 , +0.75 , +1.00 ]);

ylim([-5,21]);
ylabel('$y$','Interpreter','latex');

set(gca,'TickLabelInterpreter','latex');
set(gca,'XMinorTick','on','YMinorTick','on');





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


