function [yLast,gLast,tLast,TestFlag,TestFlagValue,n10,n90,x1,x2,alarm,tt,ft]=AggKineticsV15PhaseExplor(DataArray)
%
% ft is the quantity that is being tested. change it at the end of the
% function
%

    
% change random numbers stream
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

% get data from the DataArray

klpf = DataArray(1); 
klmf = DataArray(2);
klpm = DataArray(3);
klmm = DataArray(4);
kap = DataArray(5);
kam = DataArray(6);
kn = DataArray(7);
cnmRatio = DataArray(8);        % cnmRatio = (cnm equivalent)/cn = cnm*p0/cn
kbf = DataArray(9);
kbm = DataArray(10);
y00 = DataArray(11);
p00 = DataArray(12);
rs = DataArray(13);
rfm = DataArray(14);
n = DataArray(15);
m = DataArray(16);

% Normalization factor
% Don't forget to adjust real time accordingly

% Convert to Gillespie and Normalize parameters
%(klpf,klmf,klpm,klmm,kap,kam,kn,knm,kbf,kbm,y0,p0,rs,n,m,v0)
v0 = 1e15;          % in nm^3
na = 6.022e23;      % Avogadro Num
v1 = v0*1e-24;      % Convert volume from (nm)^3 to liters
v = v1*na;          % facilitate conversion of # of mols to # of particles
alarm  = 0;

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

% automating the scaling
vScale = 1;
% minN = min(y0,p0);
% if minN < 5e5
%     vScale = minN/1e6;      % increase the smaller number of particles to 1e6
%     alarm = 1;              % Raise scale down alarm
% elseif minN > 2e6
%     vScale = 1e6/minN;      % reduce the number of particles if the minN is larger than 2e6
%     alarm = 2;              % Raise scale up alarm
% end



% Scaling
clpf = clpf*vScale;         % a2 = clpf*y(1)*y(2:(n-1));
% clmf Stay the same        % a4 = 2*clmf*y(2:n);
clpm = clpm*vScale;         % a6 = clpm*y(1)*g(1:(n-1));
% clmm = Stay the same      % a8 = clmm*g(2:n);
cap = cap*vScale;           % a5 = cap*g0*y;
% cam Stay the same         % a9 = cam*g(1:n)
% cn : update at the calculation        % a1 = cn*y(1)*(vScale*y(1)-1)/2;  
% cbf Stay the same                     % a3 = cbf*y(4:n).*(i-3);     
% cbm Stay the same                     % a7 = cbm*g(3:n).*(i-2)    
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
% inagg = zeros(1,m);     % Collect data on sequestered aggregates
y = zeros(1,n);         % The number of free aggregates of different length
g = zeros(1,n);         % The number of mambrane aggregates of different lengt
dyin = 0;               % initialize parameters for the deterministic process
dg0in = 0;              % initialize parameters for the deterministic process
aa = zeros(1,12);       % Collect the propensity of each region
% dg0 = 0;                % Collect the fractions of the g0
% dy1 = 0;                % Collect the fractions of the y1
crfmArray = crfm./(1:n);% out flow aggregates rate is adjusted by aggregate length


% % for testing
% ynt = zeros(1,m);
% ynt(1) = y0;


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


% TestFlag = 0;       % Counts the number of time steps that the test was positive
if y00<0 || p00<0
    PositiveFlag = 0;
else
    PositiveFlag = 1;
end

if PositiveFlag
    i = 2;
    isteps = 1;
    while i<=m && isteps<1e9
        % Add time step
        if a0 ~= 0
            dt = log(1/rand)/a0;
        else
            dt = 1;
        end

        % Get the reaction type number                           
        r = rand;           % Random number (0,1)
        rn = 1;             % the reaction number in a (1 to 12)
        s = aa(rn)/a0;      % accomulating the propensities
        while s < r
            rn = rn + 1;
            s = s + aa(rn)/a0;  
        end
        % find the reaction within the reaction type
        % agn1 is the possition in the y or g arrays
        agn1 = rst(rn);        % the first position in the correct region of the propencity array
        if rn~=1 && rn~=10     % for rn=1 and rn=10 we only have one value in the propencity array                        
            s = s - aa(rn)/a0 + a(rn,agn1)/a0;      % return value to begining of reaction
            while s < r   
                agn1 =agn1 + 1;
                s = s + a(rn,agn1)/a0; 
            end
        end

        % The deterministic flow
        % Adjust y(1) and g0 
        if crfm>0
            dy1 = (crfp/crfm*(1-exp(-crfm*dt)) + y(1)*exp(-crfm*dt)) + dyin - y(1);   % for y(1)
        else
            dy1 = 0;
        end
        dy1Int = fix(dy1);                                                      % Round toward zero
        dyin = dy1 - dy1Int;                                                    % Collecting fractions of particles
        y(1) = y(1) + dy1Int;                                                  % Changing the number of particles 
        y(1) = y(1)*(y(1)>0);                                                % No negative monomer numbers

        % Membrane protein production
        if crs>0
            dg0 = (crp/crs*(1-exp(-crs*dt)) + g0*exp(-crs*dt)) + dg0in - g0;  % for g0
        else
            dg0 = 0;
        end
        dg0Int = fix(dg0);                                              % Round toward zero
        dg0in = dg0 - dg0Int;                                           % Collecting fractions of particles
        g0 = g0 + dg0Int;                                               % Changing the number of particles 
        g0 = g0*(g0>0);                                                 % No negative monomer numbers


        % Calculate the changes acording to the reaction that took place   
        switch rn
            case 1
                % Calculate y and g
                if y(1)>2
                    y(1) = y(1)-2;
                    y(2) = y(2)+1;
                    F(i) = F(i-1)+ 2;
                    % reaction propensities that need update
                    % a1 = cn*y(1)*(y(1)-1)/2;  Updated at the end
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a4 = 2*clmf*y(2:n);
                    a(4,2) = clmf*y(2);
                    aa(4) = sum(a(4,:));
                    % a5 = cap*g0*y;    Updated at the end
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a10 = cnm*g0*y(1)*(y(1)-1)/2      Updated at the end
                else
                    F(i) = F(i-1);
                end
                M(i) = M(i-1); 
            case 2
                if y(1)>0 && y(agn1)>0
                    % agn2 is the possition in the y or g arrays
                    agn2 = agn1+1;
                    % Calculate y and g
                    y(agn1) = y(agn1)-1;
                    y(1) = y(1)-1;
                    y(agn2) = y(agn2)+1;
                    F(i) = F(i-1)+ 1;
                    % Calculate changes
                    if agn1>3
                        % a3 = cbf*y(4:n).*(i-3);
                        a(3,agn1) = cbf*y(agn1)*(agn1-3);
                        a(3,agn2) = cbf*y(agn2)*(agn2-3);
                        aa(3) = sum(a(3,:));
                    elseif agn2>3
                        % a3 = cbf*y(4:n).*(i-3);
                        a(3,agn2) = cbf*y(agn2)*(agn2-3);
                        aa(3) = sum(a(3,:));
                    end
                    % reaction propensities that need update
                    % a1 = cn*y(1)*(y(1)-1)/2;     Updated at the end
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a4 = 2*clmf*y(2:n);
                    a(4,agn1) = (1+(agn1>2))*clmf*y(agn1);
                    a(4,agn2) = 2*clmf*y(agn2);
                    aa(4) = sum(a(4,:));
                    % a5 = cap*g0*y;    Updated at the end
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a10 = cnm*g0*y(1)*(y(1)-1)/2      Updated at the end
                    % a12 = crfmArray(2:n)*y(2:n);
                    a(12,agn1) = crfmArray(agn1)*y(agn1);
                    a(12,agn2) = crfmArray(agn1)*y(agn2);
                    aa(12) = sum(a(12,:));
                else
                    F(i) = F(i-1);
                end
                M(i) = M(i-1);  
            case 3
                if y(agn1)>0
                    % agn2,3 are the possition in the y or g arrays
                    agn2 = randi(agn1-3)+1;       % The aggregate can break at any bond, randomly
                    agn3 = agn1-agn2;
                    % Calculate y ang g 
                    y(agn1) = y(agn1)-1;               
                    y(agn2) = y(agn2)+1;        % y(ii)-> y(jj)+y(11-jj)
                    y(agn3) = y(agn3)+1;
                    % Calculate changes
                    if agn2~=agn3
                        % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                        % a4 = 2*clmf*y(2:n);
                        a(4,agn3) = (1+(agn3>2))*clmf*y(agn3);
                        aa(4) = sum(a(4,:));
                        % a5 = cap*g0*y;    Updated at the end
                        if agn3>3
                            % a3 = cbf*y(4:n).*yf;
                            a(3,agn3) = cbf*y(agn3)*(agn3-3);
                            aa(3) = sum(a(3,:));
                        end
                    end 
                    if agn2>3
                        % a3 = cbf*y(4:n).*yf;
                        a(3,agn2) = cbf*y(agn2)*(agn2-3);
                        aa(3) = sum(a(3,:));
                    end          
                    % reaction propensities that need update
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a3 = cbf*y(4:n).*yf; 
                    a(3,agn1) = cbf*y(agn1)*(agn1-3);
                    aa(3) = sum(a(3,:));
                    % a4 = 2*clmf*y(2:n);
                    a(4,agn1) = 2*clmf*y(agn1);
                    a(4,agn2) = (1+(agn2>2))*clmf*y(agn2);
                    aa(4) = sum(a(4,:));
                    % a5 = cap*g0*y;    Updated at the end
                    % a12 = crfmArray(2:n)*y(2:n);
                    a(12,agn1) = crfmArray(agn1)*y(agn1);
                    a(12,agn2) = crfmArray(agn2)*y(agn2);
                    a(12,agn3) = crfmArray(agn3)*y(agn3);
                    aa(12) = sum(a(12,:));
                end
                % F(i) does not change
                F(i) = F(i-1);
                M(i) = M(i-1);       
            case 4
                if y(agn1)>0
                    % agn2 is the possition in the y or g arrays
                    agn2 = agn1 - 1;
                    % Calculate y and g
                    y(agn1) = y(agn1)- 1;   
                    y(1) = y(1)+ 1;
                    y(agn2) = y(agn2)+ 1; 
                    if agn1>2
                        F(i) = F(i-1)- 1;
                    else
                        F(i) = F(i-1)- 2;
                    end
                    % Calculate changes
                    if agn1>3
                        % a3 = cbf*y(4:n).*yf;
                        a(3,agn1) = cbf*y(agn1)*(agn1-3);
                        aa(3) = sum(a(3,:));
                    end
                    if agn2>3
                        % a3 = cbf*y(4:n).*yf;
                        a(3,agn2) = cbf*y(agn2)*(agn2-3);
                        aa(3) = sum(a(3,:));
                    end
                    if agn2>1
                        % a4 = 2*clmf*y(2:n);
                        a(4,agn2) = (1+(agn2>2))*clmf*y(agn2);
                        aa(4) = sum(a(4,:));
                        % a5 = cap*g0*y;    Updated at the end
                    end
                    % reaction propensities that need update 
                    % a4 = 2*clmf*y(2:n);
                    a(4,agn1) = (1+(agn1>2))*clmf*y(agn1);
                    aa(4) = sum(a(4,:));
                    % a1 = cn*y(1)*(y(1)-1)/2;   Updated at the end
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a5 = cap*g0*y;    Updated at the end
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a10 = cnm*g0*y(1)*(y(1)-1)/2      Updated at the end
                    % a12 = crfmArray(2:n)*y(2:n);
                    a(12,agn1) = crfmArray(agn1)*y(agn1);
                    if agn2>1
                        a(12,agn2) = crfmArray(agn2)*y(agn2);
                    end
                    aa(12) = sum(a(12,:));
                else
                    F(i) = F(i-1);
                end
                M(i) = M(i-1);    
            case 5
                if y(agn1)>0 && g0>0
                    %Adjust y and g
                    y(agn1) = y(agn1)- 1;
                    g0 = g0 - 1;
                    g(agn1) = g(agn1)+ 1;              
                    F(i) = F(i-1)-agn1;  % agn1 should be always >1
                    M(i) = M(i-1)+agn1; 
                    % Calculate Changes         
                    if agn1>3
                        % a3 = cbf*y(4:n).*yf;
                        a(3,agn1) = cbf*y(agn1)*(agn1-3);
                        aa(3) = sum(a(3,:));
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn1)= cbm*g(agn1)*(agn1-2);
                        aa(7) = sum(a(7,:));
                    end
                    % reaction propensities that need update
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a4 = 2*clmf*y(2:n);
                    a(4,agn1) = (1+(agn1>2))*clmf*y(agn1);
                    aa(4) = sum(a(4,:));
                    % a5 = cap*g0*y;    Updated at the end
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a8 = clmm*g(2:n);
                    a(8,agn1)= clmm*g(agn1);
                    aa(8) = sum(a(8,:));
                    % a9 = cam*g(2:n);
                    a(9,agn1)= cam*g(agn1);
                    aa(9) = sum(a(9,:));
                    % a10 = cnm*g0*y(1)*(y(1)-1)/2      Updated at the end
                    % a11 = crs*g(1:n)
                    a(11,agn1)= crs*g(agn1);
                    aa(11) = sum(a(11,:));
                    % a12 = crfmArray(2:n)*y(2:n);
                    a(12,agn1) = crfmArray(agn1)*y(agn1);
                    aa(12) = sum(a(12,:));
                else
                    F(i) = F(i-1);  
                    M(i) = M(i-1);
                end  
            case 6
                if g(agn1)>0 && y(1)>0
                    % agn2 is the possition in the y or g arrays
                    agn2 = agn1 + 1;
                    % Calculate y and g
                    g(agn1) = g(agn1)- 1;
                    y(1) = y(1)- 1;
                    g(agn2) = g(agn2)+ 1;
                    M(i) = M(i-1)+ 1;   
                    % Calculate changes
                    if agn1>3
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn1)= cbm*g(agn1)*(agn1-2);
                        a(7,agn2)= cbm*g(agn2)*(agn2-2);
                        aa(7) = sum(a(7,:));
                    elseif agn2>3
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn2)= cbm*g(agn2)*(agn2-2);
                        aa(7) = sum(a(7,:));
                    end
                    if agn1>1
                        % a8 = clmm*g(2:n);
                        a(8,agn1)= clmm*g(agn1);
                        aa(8) = sum(a(8,:)); 
                    end
                    % reaction propensities that need update
                    % a1 = cn*y(1)*(y(1)-1)/2;     Updated at the end
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a8 = clmm*g(2:n);  agn2>1
                    a(8,agn2)= clmm*g(agn2);
                    aa(8) = sum(a(8,:));
                    % a9 = cam*g(2:n);
                    a(9,agn1)= cam*g(agn1);
                    a(9,agn2)= cam*g(agn2);
                    aa(9) = sum(a(9,:));
                    % a10 = cnm*g0*y(1)*(y(1)-1)/2      Updated at the end
                    % a11 = crs*g(1:n)
                    a(11,agn1)= crs*g(agn1);
                    a(11,agn2)= crs*g(agn2);
                    aa(11) = sum(a(11,:));
                else
                    M(i) = M(i-1);
                end
                F(i) = F(i-1);
            case 7
                if g(agn1)>0
                    % agn2,3 are the possition in the y or g arrays
                    agn2 = randi(agn1-3)+1;       % The aggregate can break at any bond, randomly (g length)
                    agn3 = agn1-agn2;           % ang3 relate to the y length
                    % Calculate y and g
                    g(agn1) = g(agn1)- 1; 
                    g(agn2) = g(agn2)+ 1;
                    y(agn3) = y(agn3)+ 1;  
                    F(i) = F(i-1)+agn3;
                    M(i) = M(i-1)-agn3; 
                    % Calculate changes
                    if agn3>3
                        % a3 = cbf*y(4:n).*yf;
                        a(3,agn3) = cbf*y(agn3)*(agn3-3);
                        aa(3) = sum(a(3,:));
                    end
                    if agn3>1
                        % a4 = 2*clmf*y(2:n);
                        a(4,agn3) = (1+(agn3>2))*clmf*y(agn3);
                        aa(4) = sum(a(4,:));
                    end
                    if agn2>1
                        % a8 = clmm*g(2:n);
                        a(8,agn2)= clmm*g(agn2);
                        aa(8) = sum(a(8,:));
                    end
                    if agn2>3
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn2)= cbm*g(agn2)*(agn2-2);
                        aa(7) = sum(a(7,:));
                    end  
                    % reaction propensities that need update
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end  
                    % a5 = cap*g0*y;    Updated at the end
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a7 = cbm*g(4:n).*(i-2)
                    a(7,agn1)= cbm*g(agn1)*(agn1-2);
                    aa(7) = sum(a(7,:));
                    % a8 = clmm*g(2:n);
                    % a8 = clmm*g(2:n);
                    a(8,agn1)= clmm*g(agn1);
                    aa(8) = sum(a(8,:));
                    % a9 = cam*g(2:n);
                    a(9,agn1)= cam*g(agn1);
                    a(9,agn2)= cam*g(agn2);
                    aa(9) = sum(a(9,:));
                    % a11 = crs*g(1:n)
                    a(11,agn1)= crs*g(agn1);
                    a(11,agn2)= crs*g(agn2);
                    aa(11) = sum(a(11,:));
                    % a12 = crfmArray(2:n)*y(2:n);
                    if agn3>1
                        a(12,agn3) = crfmArray(agn3)*y(agn3);
                    end
                    aa(12) = sum(a(12,:));
                else
                    F(i) = F(i-1);
                    M(i) = M(i-1);
                end 
            case 8
                if g(agn1)>0
                    % agn2 is the possition in the y or g arrays
                    agn2 = agn1 - 1; 
                    % Calculate y and g
                    if agn1>2
                        g(agn1) = g(agn1)- 1;   
                        g(agn2) = g(agn2)+ 1; 
                        y(1) = y(1)+ 1;
                        M(i) = M(i-1)- 1;   
                    else
                        % if agn1=2 it falls apart to 2 seperate monomers
                        g(agn1) = g(agn1)- 1;   
                        g0 = g0+ 1; 
                        y(1) = y(1)+ 2;
                        M(i) = M(i-1)- 2; 
                    end
                    % Calculate changes
                    if agn1>3
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn1)= cbm*g(agn1)*(agn1-2);
                        aa(7) = sum(a(7,:));
                    end
                    if agn2>3
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn2)= cbm*g(agn2)*(agn2-2);
                        aa(7) = sum(a(7,:));
                    end
                    if agn2>1
                        % a8 = clmm*g(2:n);
                        a(8,agn2)= clmm*g(agn2);
                        aa(8) = sum(a(8,:));
                    end
                    % reaction propensities that need update
                    % a1 = cn*y(1)*(y(1)-1)/2;     Updated at the end
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end    
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end  
                    % a8 = clmm*g(2:n);
                    a(8,agn1)= clmm*g(agn1);
                    aa(8) = sum(a(8,:));
                    % a9 = cam*g(2:n);
                    a(9,agn1)= cam*g(agn1);
                    if agn2>1
                        a(9,agn2)= cam*g(agn2);
                    end
                    aa(9) = sum(a(9,:));
                    % a10 = cnm*g0*y(1)*(y(1)-1)/2      Updated at the end
                    % a11 = crs*g(1:n)
                    a(11,agn1)= crs*g(agn1);
                    if agn2>1
                        a(11,agn2)= crs*g(agn2);
                    end
                    aa(11) = sum(a(11,:));
                else
                    M(i) = M(i-1);
                end
                F(i) = F(i-1);
            case 9
                if  g(agn1)>0
                    % Calculate y and g
                    y(agn1) = y(agn1)+ 1;
                    g0 = g0 + 1;
                    g(agn1) = g(agn1)- 1;
                    if agn1>1
                        F(i) = F(i-1)+ agn1;
                    else
                        F(i) = F(i-1);
                    end
                    M(i) = M(i-1)- agn1;
                    % Calculate changes
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a4 = 2*clmf*y(2:n);
                    a(4,agn1) = (1+(agn1>2))*clmf*y(agn1);
                    aa(4) = sum(a(4,:));
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a8 = clmm*g(2:n);
                    a(8,agn1)= clmm*g(agn1);
                    aa(8) = sum(a(8,:));
                    if agn1>3
                        % a3 = cbf*y(4:n).*yf;
                        a(3,agn1) = cbf*y(agn1)*(agn1-3);
                        aa(3) = sum(a(3,:));
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn1)= cbm*g(agn1)*(agn1-2);
                        aa(7) = sum(a(7,:));
                    end         
                    % reaction propensities that need update
                    % a1 = cn*y(1)*(y(1)-1)/2;     Updated at the end
                    % a5 = cap*g0*y;    Updated at the end
                    % a9 = cam*g(2:n);
                    a(9,agn1)= cam*g(agn1);
                    aa(9) = sum(a(9,:));
                    % a10 = cnm*g0*y(1)*(y(1)-1)/2      Updated at the end
                    % a11 = crs*g(1:n)
                    a(11,agn1)= crs*g(agn1);
                    aa(11) = sum(a(11,:));
                    % a12 = crfmArray(2:n)*y(2:n);
                    a(12,agn1) = crfmArray(agn1)*y(agn1);
                    aa(12) = sum(a(12,:));
                else
                    F(i) = F(i-1);
                    M(i) = M(i-1);
                end
            case 10
                if y(1)>0 && g0>0 && cnm>0
                    %Adjust y and g
                    y(1) = y(1) - 2;
                    g0 = g0 - 1;
                    g(2) = g(2)+ 1;     
                    M(i) = M(i-1) + 2;   
                    % Calculate changes
                    % reaction propensities that need update
                    % a1 = cn*y(1)*(y(1)-1)/2;      Updated at the end
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a5 = cap*g0*y;    Updated at the end
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a8 = clmm*g(2:n);
                    a(8,2)= clmm*g(2);
                    aa(8) = sum(a(8,:));
                    % a9 = cam*g(2:n);
                    a(9,2)= cam*g(2);
                    aa(9) = sum(a(9,:));
                    % a11 = crs*g(1:n)
                    a(11,2)= crs*g(2);
                    aa(11) = sum(a(11,:));
                else
                    M(i) = M(i-1);
                end
                F(i) = F(i-1);      % monomers are not counted in F
            case 11
                if g(agn1)>0
                    g(agn1) = g(agn1) - 1;
                    M(i-1) = M(i-1) - agn1;
                    % Do not do anything to F since i did not change
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    if agn1>3
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn1)= cbm*g(agn1)*(agn1-2);
                        aa(7) = sum(a(7,:));
                    end  
                    % a8 = clmm*g(2:n);
                    a(8,agn1)= clmm*g(agn1);
                    aa(8) = sum(a(8,:));
                    % a9 = cam*g(2:n);
                    a(9,agn1)= cam*g(agn1);
                    aa(9) = sum(a(9,:));
                    % a11 = crs*g(1:n)
                    a(11,agn1)= crs*g(agn1);
                    aa(11) = sum(a(11,:));
                end
            case 12
                if y(agn1)>0
                    y(agn1) = y(agn1) - 1;
                    F(i-1) = F(i-1) - agn1;
                    % Do not do anything to M since i did not change
                    % Calculate changes
                    % reaction propensities that need update
                    % a1 = cn*y(1)*(y(1)-1)/2;      Updated at the end
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end
                    if agn1>3
                        % a3 = cbf*y(4:n).*(i-3);
                        a(3,agn1) = cbf*y(agn1)*(agn1-3);
                        aa(3) = sum(a(3,:));
                    end
                    % a4 = 2*clmf*y(2:n);
                    a(4,agn1) = (1+(agn1>2))*clmf*y(agn1);
                    aa(4) = sum(a(4,:));
                    % a5 = cap*g0*y;    Updated at the end
                    % a12 = crfmArray(2:n)*y(2:n);
                    a(12,agn1) = crfmArray(agn1)*y(agn1);
                    aa(12) = sum(a(12,:));  
                end    
            otherwise
                % disp('Problem with the reaction number')   
        end

        if a(11,1)~=0
            disp('a(11,1)~=0')
            break
        end

        % Always update 
        y0t(i) = y(1);
        g0t(i) = g0;
        % reaction propensities that need update
        % a1 = cn*y(1)*(y(1)-1)/2;  
        a(1,1) = cn*y(1)*(vScale*y(1)-1)/2;
        aa(1) = a(1,1);
        % a2 = clpf*y(1)*y(2:(n-1));
        a(2,rst(2):re(2))= clpf*y(1)*y(2:(n-1));
        aa(2) = sum(a(2,:));
        % a5 = cap*g0*y;
        a(5,rst(5):re(5))= cap*g0*y(2:end);
        aa(5) = sum(a(5,:));
        % a6 = clpm*y(1)*g(1:(n-1));
        a(6,rst(6):re(6))= clpm*y(1)*g(2:(n-1));
        aa(6) = sum(a(6,:));
        % a10 = cnm*g0*y(1)*(y(1)-1)/2
        a(10,1) = cnm*g0*y(1)*(vScale*y(1)-1)/2;
        aa(10) = a(10,1); 
        % Sum propensities
        a0 = sum(aa);

        % adjust time steps
        if rn<11
            % Aggregation step
            t(i) = t(i-1)+ dt; 
            i = i + 1;
            isteps = isteps + 1;
        else
            % clearance step
            % do not increase the i counter.
            t(i-1) = t(i-1)+ dt;
            isteps = isteps + 1;
        end
    end
    % test results
    % Test how many of the last 500 steps had segnificant aggregation    
    %  TestVal = 100*F(end-500:end)/(F(end-500:end)+y0t(end-500:end));   
    [TestFlag,TestFlagValue,n10,n90,x1,x2] = TestABetaValue(y0t,t);
%     [TestFlag,TestFlagValue,n10,n90,x1,x2] = TestAggValue(F,t);
    %  I take out the last term since in high flow rates it can cause
    %  problems
else
    TestFlag = 1;                   % Test fails -> TestFlag=1
    TestFlagValue = 1;
end



% Real world time: treal = t(m)*v/3600;
t = t/3600/24;      % Covert Gillespie time to real time 
% ft = F*vScale/v;
ft = y0t*vScale/v;
% convert concentration
% F = F*vScale/v;
% M = M*vScale/v;
% y0t = y0t*vScale/v;
% g0t = g0t*vScale/v;

% Reduce array size
tt = [t(1:9) t(10:10:end)];
ft = [ft(1:9) ft(10:10:end)];
n90 = find(tt>t(n90),1);
n10 = find(tt>t(n10),1);

yLast = y(end);
gLast = g(end);
tLast = t(end);
end