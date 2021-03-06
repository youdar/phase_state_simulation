
% General initialization code for the aggregation functions
% Use by eval('AggKineticsV16Input')

% Free aggregates
% a1 = cn*y(1)*(y(1)-1)/2;
% a2 = clpf*y(1)*y(2:(n-1));
% a3 = cbf*y(4:n).*(i-3);
% a4 = 2*clmf*y(2:n);
% Membrane aggregates
% a5 = cap*g0*y;
% a6 = clpm*y(1)*g(2:(n-1));
% a7 = cbm*g(4:n).*(i-2)
% a8 = clmm*g(2:n);
% a9 = cam*g(2:n);
% a10 = cnm*g0*y(1)*(y(1)-1)/2
% a11 = crs*g(2:n)
% a12 = crfmArray(2:n)*y(2:n);    crfmArray= crfm./(2:n)
% production and sequestration of membrane and free monomers is handled
% deterministically

% klpf = 1.5e3 ;
% klmf = 3e-8;
% klpm = 1.5e3 ;
% klmm = 3e-8;
% kap = 100;
% kam = 5e-7;
% kn = 9.5e-5;
% cnmRatio = 1;
% kbf = 1.5e-9;
% kbm = 1.5e-9;    
% y00 = 2e-10; 
% p00 = 2e-9;  
% rs = 1e-8;
% rfm =1e-8; 
% n = 500;
% m = 5e4;

% klpf = 3e8;
% klmf = 2.45e-2;
% klpm = 3e8;
% klmm = 2.45e-2;
% kap = 0.5e5;
% kam = 2.5e-4;
% kn = 2e-5;
% cnmRatio = 1;
% kbf = 1e-6;
% kbm = 1e-6;    
% y00 = 2e-10; 
% p00 = 2e-9;  
% rs = 2.5e-5;
% rfm = 2e-5; 
% n = 500;
% m = 1e6;


% Use the basic data from figure 4.B and the Phy. Rev. E paper
klpf = 1.5e3 ;
klmf = 3e-8;
klpm = 1.5e3 ;
klmm = 3e-8;
kap = 100;
kam = 5e-7;
kn = 9.5e-5;
cnmRatio = 1;
kbf = 1.5e-9;
kbm = 1.5e-9;    
y00 = 2e-9; 
p00 = 2e-10;  
rs = 2e-8;
rfm =2e-8; 
n = 500;
m = 1.5e6;
