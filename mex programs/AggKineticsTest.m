% function AggKineticsV7(klpf,klmf,klpm,klmm,kap,kam,kn,knm,kbf,kbm,y00,p00,rs,rfm,n,m,v0)
% Using a combination of Gillespie algorithm and deterministc solutions
% to investigate protein aggregation kinetics
%
% Plots will be saved in current directory
%
% Some General Data:
% Prion cycle ~60 min (Prion biology & diseas. p 493)
% Prion particle per 1 micro meter of membrane table (Prion biology & diseas. p 490)   
%     Startum oriens - axon axolemma  ~0.95 (0.00095 per 1 nm of membrane)
% Ultrathin Cyrosection, brain slice thinkness ~ 100 nm 
% We get ~95,000 PrPc/nm^2
% Area of neuron cell membrane ~ 250,000 nm^2
%
% v0[qubic micro meters]
% k's reaction rate constants in [1/(M*s)], [1/s]...
% y00,p00 in [M]
%
% klpf k+ for elongation (monomer attachment) for free aggregates
% klmf k- for monomer dettachment, for free aggregates
% klpm k+ for elongation (monomer attachment) for membrane bounnd aggregates
% klmm k- for monomer dettachment, for membrane bounnd aggregates 
% kap/kam attachemnt and dettachment rates, of aggregates to te membrane
% kn free aggregates nucleation rate
% knm mambrane aggregation rate
% kbf free aggregates brakage rate
% kbm mambrane aggregates brakage rate
% y00 initial concentration of free monomers (A-beta)
% p00 initial concentration of membrane bound proteins (Prion)
% rs sequestration rate of membrane bound proteins and aggregates
% rfp new, free monomers supply rate
% rp membrane bound proteins production rate
% bn burst size, the numbers of new proteins that enters to the system at
%                              each production step
% n the max length of aggregates
% m the number of iterations we want to do 
% v0 the system volume. used to convert concentration and reaction rates to
%                             stochastic form




clc         % Clear work space
close       % Close all figures
% (klpf,klmf,klpm,klmm,kap,kam,kn,knm,kbf,kbm,y00,p00,rs,rfm,n,m,v0)
% (3e3,2e-8,3e3,2e-8,4e1,2e-8,2e-4,2e-5,2e-8,2e-10,2e-9,10e-9,5e-9,5e-6,800,9e6,1e6);
% (3e4,2e-6,3e4,2e-6,1e2,2e-7,2e-5,2e7,2e-8,2e-10,1e-10,2e-9,1e-5,1e-5,500,5e6,1e6);
clear       % Clear all variables
klpf = 6e4;
klmf = 2e-6;
klpm = 6e4;
klmm = 2e-6;
kap = 1e2;
kam = 2e-7;
kn = 6e-5;
knm = 2e7;
kbf = 2e-8;
kbm = 2e-10;
y00 = 1e-10;
p00 = 2e-9;
rs = 1e-5;
rfm = 1e-7;
n = 500;
m = 1e6;
v0 = 1e6;




% Convert and Normalize Gillespie parameters
%(klpf,klmf,klpm,klmm,kap,kam,kn,knm,kbf,kbm,y0,p0,rs,n,m,v0)
na = 6.022e23;      % Avogadro Num
v1 = v0*1e-15;       % Convert volume to liters
v = v1*na;           % facilitate conversion of # of mols to # of particles

clpf = 1;
clmf = v*klmf/klpf;
clpm = klpm/klpf;
clmm = v*klmm/klpf;
cap = kap/klpf;
cam = v*kam/klpf;
cn = 2*kn/klpf;         % free nucleation
cnm = 2*knm/klpf/v;     % nucleation on the membrane
cbf = v*kbf/klpf;       % Free aggregates breaking rate
cbm = v*kbm/klpf;       % Membrane aggregates breaking rate
y0 = round(v*y00);      % Convert concentration to number of particles in v0
p0 = round(v*p00);      % Convert concentration to number of particles in v0
rp = rs*p00;
crs = v*rs/klpf;
crp = v^2*rp/klpf;   % crp = crs*p0
rfp = rfm*y00;          % flow supply rate
crfm = v*rfm/klpf;
crfp = crfm*y0; % flow supply of free monomers (crfp = crfm*y0)v^2*rfm*y00/klpf = v*rfm*y0/klpf


% Display some numbers:
disp ('~~~~~~~~~~~~~~~~~~~~~~~')
disp (['klpf = ' num2str(klpf) '               klmf = ' num2str(klmf)])
disp (['klpm = ' num2str(klpm) '               klmm = ' num2str(klmm)])
disp (['kap = ' num2str(kap) '                kam = ' num2str(kam)])
disp (['kn = ' num2str(kn)])
disp (['knm = ' num2str(knm)])
disp (['kbf = ' num2str(kbf)])
disp (['kbm = ' num2str(kbm)])
disp (['p00 = ' num2str(p00)])
disp (['y00 = ' num2str(y00)])
disp (['rs = ' num2str(rs) '                  rp = ' num2str(rp)])
disp (['rfm = ' num2str(rfm) '                 rfp = ' num2str(rfp)])
disp ('=========================')
disp (['clpf = ' num2str(clpf) '               clmf = ' num2str(clmf)])
disp (['clpm = ' num2str(clpm) '               clmm = ' num2str(clmm)])
disp (['cap = ' num2str(cap) '                cam = ' num2str(cam)])
disp (['cn = ' num2str(cn)])
disp (['cnm = ' num2str(cnm)])
disp (['cbf = ' num2str(cbf)])
disp (['cbm = ' num2str(cbm)])
disp (['p0 = ' num2str(p0)])
disp (['y0 = ' num2str(y0)])
disp (['crs = ' num2str(crs) '                crp = ' num2str(crp)])
disp (['crfm = ' num2str(crfm) '               crfp = ' num2str(crfp)])
disp ('~~~~~~~~~~~~~~~~~~~~~~~')


% Initialize parameters

t = zeros(1,m);         % collect the time steps
F = zeros(1,m);         % Collect the amount of monomers in the free aggregates at each time step
M = zeros(1,m);         % Collect the amount of monomers in the membrane aggregates at each time step
y0t = zeros(1,m);       % the number of free monomers
g0t = zeros(1,m);       % the number of membrane single proteins
inagg = zeros(1,m);     % Collect data on sequestered aggregates
y = zeros(1,n);         % The number of free aggregates of different length
g = zeros(1,n);         % The number of mambrane aggregates of different lengt
dyin = 0;               % initialize parameters for the deterministic process
dg0in = 0;              % initialize parameters for the deterministic process
aa = zeros(1,10);       % Collect the propensity of each region
dfy = zeros(1,n);       % Collect the fractions of the specific y
dfg = zeros(1,n);       % Collect the fractions of the specific g
dfg0 = 0;               % Collect the fractions of the g0



% nc = 2;               % Critical nucleation length - for membrane attachment
y(1) = y0;              % initial amount of monomers
yn = y(1);              % track the total number of free monomers and aggregates
g0 = p0;                % initial amount of monomers 
gn = g0;                 % track the total number of membrane aggregates
y0t(1)=y(1);            % tracking free monomers as a function of time
g0t(1)=g0;              % tracking membrane bound monomers as a function of time
% rst = [1 2 4 2 2 1 3 2 1 1];           % the start i in each reaction region
% re = [1 (n-1) n n n (n-1) n n n 1];   % the end i in each reaction region
rst = [0 1 3 1 1 0 2 1 0 0];          % the start i in each reaction region
re = [0 (n-2) (n-1) (n-1) (n-1) (n-2) (n-1) (n-1) (n-1) 0];   % the end i in each reaction region


% initial values for the different reaction propensities
% Do not include in/out-flow reactions in the Gillespie process  
% a(1),a(10) are length 1
a = zeros(10,800);
      
% None zero initial values
a(1,1) = cn*y(1)*(y(1)-1)/2;
a(10,1) = cnm*g0*y(1)*(y(1)-1)/2;



% Sum all the none zero values
aa(1) = a(1,1);
aa(10) = a(10,1);
a0 = sum(aa);

disp ('xxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
disp (['Free nucleation weight = ' num2str(aa(1))])
disp (['On memebrane nucleation weight = ' num2str(aa(10))])
disp ('xxxxxxxxxxxxxxxxxxxxxxxxxxxxx')

% "a" terms calculation
    % Free aggregates
    % a1 = cn*y(1)*(y(1)-1)/2;
    % a2 = clpf*y(1)*y(2:(n-1));
    % a3 = cbf*y(4:n).*yf;
    % a4 = 2*clmf*y(2:n);
    % Membrane aggregates
    % a5 = cap*g0*y;
    % a6 = clpm*y(1)*g(1:(n-1));
    % a7 = cbm*g(3:n).*(i-2)
    % a8 = clmm*g(2:n);
    % a9 = cam*g(1:n);
    % a10 = cnm*g0*y(1)*(y(1)-1)/2
    % Membrane protein sequstration - seperate process
    % Membrane aggregates sequstration - seperate process    
    % Out-flow - seperate process;
    % In-flow of free monomers and membrane single proteins is done
    % deterministically. since they are fixed rates and happen constantly
    %
    % out-flow is done by selecting one size aggregate to be cleared at
    % each time step, according to the probabilty of the 

tic;
% m number of iterations


disp(['rst[5] = ' num2str(rst(6))]);
disp(['re[5] = ' num2str(re(6))]);
%test(t,F,M,y0t,g0t,y,g,inagg,a0,a,aa,rst,re,yn,gn,clpf,clmf,clpm,clmm,cap,cam,cn,cnm,cbf,cbm,y0,p0,rp,crs,crp,rfp,crfm,crfp,m,n,dfy,dfg);
test(a,rst,re,y,g,clpm);
% AggKineticLoopTest(t,F,M,y0t,g0t,y,g,inagg,a0,a,aa,rst,re,yn,gn,clpf,clmf,clpm,clmm,cap,cam,cn,cnm,cbf,cbm,y0,p0,rp,crs,crp,rfp,crfm,crfp,m,n,dfy,dfg);

rtime = toc;
disp(['Program run time' num2str(rtime)])
% treal = t(m)*v/klpf/3600;
% disp (['treal = ' num2str(treal)])
% disp ('+++++++++++++++++++++++++++++')
% disp (['y(end) = ' num2str(y(end))])
% disp (['g(end) = ' num2str(g(end))])
% disp ('+++++++++++++++++++++++++++++')


% h=figure(1);
% plot(t,M,'--',t,F,'-',t,y0t,'-.',t,g0t,':',t,inagg,'.y');
% legend('Membrane','Free','A-Beta','PrPc','Sequestered','Location','NEO');
% title('Total momomers in aggregates');
% %xlabel('t[h]');
% xlabel('t');
% ylabel('[#]');
% ymax = max([M F y0t g0t]);
% % ymax = max([M F y0t]);
% ylim([0 ymax*1.05]);
% xlim([0 t(end)]);
% % set(gca,'XTick',0:5:treal);
% % set(gca,'XTickLabel',num2str(0:5:treal));
% % (klpf,klmf,klpm,klmm,kap,kam,kn,knm,kbf,kbm,y00,p00,rs,rfm,n,m,v0)
% str1 = sprintf('klpf=%.2g, klmf=%.2g \nklpm=%.2g, klmm=%.2g ',klpf,klmf,klpm,klmm);
% str2 = sprintf('\nkap=%.2g, kam=%.2g \nkn=%.2g, knm=%.2g',kap,kam,kn,knm);
% str3 = sprintf('\nkbf=%.2g, kbm=%.2g  \ny=%.2g, p=%.2g',kbf,kbm,y00,p00);
% str4 = sprintf('\nrs=%.2g ,rp=%.2g \nrfm=%.2g, rfp=%.2g',rs,rp,rfm,rfp);
% str5 = sprintf('\nn=%.2g, m=%.2g \nv0=%.2g',n,m,v0);
% str6 = sprintf('\ny0=%.2g, p0=%.2g',y0,p0);
% str7 = sprintf('\n\nRun time: %g \nReal Time: %g',rtime,treal);
% plottext = [str1 str2 str3 str4 str5 str6 str7];
% %gtext(plottext)
% Ypos = ymax*0.4;
% Xpos = t(end)*1.05;
% text(Xpos,Ypos,plottext)
% %filename = ['AggKin plot ' datestr(clock) '.jpg'];
% %fname = sprintf('y0=%.2g p0=%.2g kbm=%.2g klmm=%.2g kam=%.2g',y0,p0,kbm,klmm,kam);
% fname = sprintf('y00=%.2g kap=%.2g rs=%.2g rfm=%.2g',y00,kap,rs,rfm);
% filename = ['AggKin ' fname '.jpg'];
% filename = strrep(filename,':','-');
% saveas(h,filename);

% figure(2)
% plot(y(2:end))
% % xlabel('y');
% % ylabel('[#]');
% title('Free aggregates distribution');
% 
% figure(3)
% plot(g)
% % xlabel('a');
% % ylabel('[#]');
% title('membrane aggregates distribution');
% 
% figure(4)
% plot(inagg)
% title('Sequestered aggregates');
