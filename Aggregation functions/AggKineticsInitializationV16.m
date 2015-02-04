% General initialization code for the aggregation functions
% Use by eval('AggKineticsInitializationV16')

% Convert to Gillespie and Normalize parameters
%(klpf,klmf,klpm,klmm,kap,kam,kn,knm,kbf,kbm,y0,p0,rs,n,m,v0)
v0 = 1e15;          % in nm^3
na = 6.022e23;      % Avogadro Num
v1 = v0*1e-24;      % Convert volume from (nm)^3 to litters 
v = v1*na;          % facilitate conversion of # of mols to # of particles

clpf = klpf/v;
clmf = klmf;
clpm = klpm/v;
clmm = klmm;
cap = kap/v;
cam = kam;
cn = 2*kn/v;        % free nucleation
cbf = kbf;          % Free aggregates breaking rate
cbm = kbm;          % Membrane aggregates breaking rate
y0 = round(v*y00);  % Convert concentration to number of particles in v0
p0 = round(v*p00);  % Convert concentration to number of particles in v0
% rp = rs*p00;
% non-Gillespie parameters
crs = rs;       % Prion Sequestration rate
crp = crs*p0;   % Production rate        
crfm = rfm;     % Free aggregates out flow rate
crfp = crfm*y0; % flow supply of free monomers (crfp = crfm*y0)v^2*rfm*y00 = v*rfm*y0


% Scalling  factor. 
vScale = 1;
clpf = clpf*vScale;         % a2 = clpf*y(1)*y(2:(n-1));
% clmf Stay the same        % a4 = 2*clmf*y(2:n);
clpm = clpm*vScale;         % a6 = clpm*y(1)*g(1:(n-1));
% clmm = Stay the same      % a8 = clmm*g(2:n);
cap = cap*vScale;           % a5 = cap*g0*y;
% cam Stay the same         % a9 = cam*g(2:n)
% cn : update at the calculation        % a1 = cn*y(1)*(vScale*y(1)-1)/2;  
% cbf Stay the same                     % a3 = cbf*y(4:n).*(i-3);     
% cbm Stay the same                     % a7 = cbm*g(4:n).*(i-2)    
y0 = round(y0/vScale);      
p0 = round(p0/vScale);     
% crs Stay the same
crp = crp/vScale;           
% crfm = Stay the same
crfp = crfp/vScale;
% t = t*v/3600;   Convert to real time

if p0~=0
    cnm = vScale*cnmRatio*cn/p0;        % Note the scaling
    % a10 = vScale*cnm*g0*y(1)*(vScale*y(1)-1)/2
else
    cnm = 0;
end


% Initialize parameters

t = zeros(1,m);         % collect the time steps
F = zeros(1,m);         % Collect the amount of monomers in the free aggregates at each time step
M = zeros(1,m);         % Collect the amount of monomers in the membrane aggregates at each time step
y0t = zeros(1,m);       % the number of free monomers
g0t = zeros(1,m);       % the number of membrane single proteins
inagg = zeros(1,m);     % Collect data on sequestered aggregates
y = zeros(1,n);         % The number of free aggregates of different length
g = zeros(1,n);         % The number of membrane aggregates of different lengt
dyin = 0;               % initialize parameters for the deterministic process
dg0in = 0;              % initialize parameters for the deterministic process
aa = zeros(1,12);       % Collect the propensity of each region
dg0 = 0;                % Collect the fractions of the g0
dy1 = 0;                % Collect the fractions of the y1
% the aggregate clearance is being modified to represent an effective
% aggregates clearance rate
crfmArray = crfm*ones(1,n);
% crfmArray(2:n) = crfmArray(2:n)/50;  
crsArray = crs*ones(1,n);
% crsArray(2:n) = crsArray(2:n)/50; 

% Critical nucleation length in this simulation is 2
y(1) = y0;              % initial amount of monomers
g0 = p0;                % initial amount of monomers 
y0t(1)=y(1);            % tracking free monomers as a function of time
g0t(1)=g0;              % tracking membrane bound monomers as a function of time
rst = [1 2 4 2 2 2 4 2 2 1 2 2];            % the start i in each reaction region
re = [1 (n-1) n n n (n-1) n n n 1 n n];     % the end i in each reaction region


% initial values for the different reaction propensities
% Do not include in/out-flow reactions in the Gillespie process  
% a(1),a(10) are length 1
a = zeros(12,n);
      
% None zero initial values
a(1,1) = cn*y(1)*(vScale*y(1)-1)/2;
a(10,1) = cnm*g0*y(1)*(vScale*y(1)-1)/2;

% Sum all the none zero values
aa(1) = a(1,1);
aa(10) = a(10,1);
a0 = sum(aa);