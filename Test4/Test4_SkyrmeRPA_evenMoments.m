% MATLAB R2018a
clear;
clc;
close all;
addpath('../functions');

Omegab = 250;

A = function_readmatd( '../Test2/skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/' , 'A.mat' );
B = function_readmatd( '../Test2/skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/' , 'B.mat' );

[ Omegas , RPA_BIS , RPA_BIV ] = function_readSkyrmeRPAeigenfreq( '../Test2/skyrme_rpa/skyrme_rpa_MODIFIED/' , 'skyrme_rpa.out' );



%% ISOSCALAR
fprintf('\nISOSCALAR\n');
F20 = function_readvecd( '../Test2/skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/' , 'F02_ISOSCALAR.vec' );
F02 = function_readvecd( '../Test2/skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/' , 'F20_ISOSCALAR.vec' );

PlotData = cell(0,1);

for kthMoment = [ 0 , 2 , 4 , 6 , 8 ]
    
    if( kthMoment == 0 )
        kernel = 'Jackson';
    else
        kernel = 'Dirichlet';
    end
    
    x = [];
    y = [];
    
    for N_it = arrayfun( @(n)2^n , [3:12] )
 
        mun = function_ChebyshevCoefficients( A , B , F20 , F02 , Omegab , N_it );
        mun = function_applyKernel( mun , kernel );

        mk_Cheb    = function_ChebyshevMoment( mun , kthMoment , Omegab );

        mk_correct = function_kthMoment( Omegas , RPA_BIS , kthMoment );

        relerr     = abs(mk_Cheb-mk_correct)/abs(mk_correct);

        fprintf('k = %1d, N_it = %5d, mk(Chebyshev) = %e, mk(correct) = %e, rel. err. = %e\n', kthMoment , N_it , mk_Cheb , mk_correct , relerr );

        x = [ x , N_it   ];
        y = [ y , relerr ];
    end
    
    PlotData{end+1} = cell( { kthMoment , x , y } );
    
    fprintf('\n');
end



figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,1,0.65,0.8]);
Legends  = cell(0);
AllMarks = {'o','d','^','s','x','.','+','*','v','>','<','p','h'};
for i = 1 : length(PlotData)
    kthMoment = PlotData{i}{1};
    N_it      = PlotData{i}{2};
    relerr    = PlotData{i}{3};
    
    Legends{i} = '$k =$ ' + string(kthMoment);
    
    semilogy( log2(N_it) , relerr , AllMarks(mod(i,length(AllMarks)))+"-" , 'LineWidth' , 2 , 'MarkerSize' , 8 ); hold on;
    
    XtickLabels = cell(0);
    for i = 1 : length(N_it)
        %XtickLabels{i} = '2^{' + string( round( log2(N_it(i)) ) ) + '}';
        XtickLabels{i} = round( N_it(i) );
    end
    xticks( log2(N_it) );
    xticklabels( XtickLabels );
end
grid on; %grid minor;

legend(Legends,'Interpreter','latex','FontSize',20);
xlabel('$N_{\mathrm{it}}$','Interpreter','latex','FontSize',30);
ylabel('$|m_k-m_k^{\mathrm{(Chebyshev)}}|/m_k$','Interpreter','latex','FontSize',30);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',20);   
    


%% ISOVECTOR
fprintf('\nISOSCALAR\n');
F20 = function_readvecd( '../Test2/skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/' , 'F02_ISOVECTOR.vec' );
F02 = function_readvecd( '../Test2/skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/' , 'F20_ISOVECTOR.vec' );

PlotData = cell(0,1);

for kthMoment = [ 0 , 2 , 4 , 6 , 8 ]
    
    if( kthMoment == 0 )
        kernel = 'Jackson';
    else
        kernel = 'Dirichlet';
    end
    
    x = [];
    y = [];
    
    for N_it = arrayfun( @(n)2^n , [3:12] )
 
        mun = function_ChebyshevCoefficients( A , B , F20 , F02 , Omegab , N_it );
        mun = function_applyKernel( mun , kernel );

        mk_Cheb    = function_ChebyshevMoment( mun , kthMoment , Omegab );

        mk_correct = function_kthMoment( Omegas , RPA_BIV , kthMoment );

        relerr     = abs(mk_Cheb-mk_correct)/abs(mk_correct);

        fprintf('k = %1d, N_it = %5d, mk(Chebyshev) = %e, mk(correct) = %e, rel. err. = %e\n', kthMoment , N_it , mk_Cheb , mk_correct , relerr );

        x = [ x , N_it   ];
        y = [ y , relerr ];
    end
    
    PlotData{end+1} = cell( { kthMoment , x , y } );
    
    fprintf('\n');
end



figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,1,0.65,0.8]);
Legends  = cell(0);
AllMarks = {'o','d','^','s','x','.','+','*','v','>','<','p','h'};
for i = 1 : length(PlotData)
    kthMoment = PlotData{i}{1};
    N_it      = PlotData{i}{2};
    relerr    = PlotData{i}{3};
    
    Legends{i} = '$k =$ ' + string(kthMoment);
    
    semilogy( log2(N_it) , relerr , AllMarks(mod(i,length(AllMarks)))+"-" , 'LineWidth' , 2 , 'MarkerSize' , 8 ); hold on;
    
    XtickLabels = cell(0);
    for i = 1 : length(N_it)
        %XtickLabels{i} = '2^{' + string( round( log2(N_it(i)) ) ) + '}';
        XtickLabels{i} = round( N_it(i) );
    end
    xticks( log2(N_it) );
    xticklabels( XtickLabels );
end
grid on; %grid minor;

legend(Legends,'Interpreter','latex','FontSize',20);
xlabel('$N_{\mathrm{it}}$','Interpreter','latex','FontSize',30);
ylabel('$|m_k-m_k^{\mathrm{(Chebyshev)}}|/m_k$','Interpreter','latex','FontSize',30);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',20);   

