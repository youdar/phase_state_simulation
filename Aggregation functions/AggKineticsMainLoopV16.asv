% Main loop used by all the aggregation functions 
% to include code use eval('AggKineticsMainLoopV16')

% m number of iterations
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
    
    
    % Calculate the changes according to the reaction that took place   
    % F is the amount of ABeta in free aggregates
    % M is the amount of ABeta in membrane aggregates
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
                % a12 = crfmArray(2:n)*y(2:n);
                a(12,2) = crfmArray(2)*y(2);
                aa(12) = sum(a(12,:));
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
                a(11,agn1)= crsArray(agn1)*g(agn1);
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
                a(11,agn1)= crsArray(agn1)*g(agn1);
                a(11,agn2)= crsArray(agn2)*g(agn2);
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
                a(11,agn1)= crsArray(agn1)*g(agn1);
                a(11,agn2)= crsArray(agn2)*g(agn2);
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
                a(11,agn1)= crsArray(agn1)*g(agn1);
                if agn2>1
                    a(11,agn2)= crsArray(agn2)*g(agn2);
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
                a(11,agn1)= crsArray(agn1)*g(agn1);
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
                a(11,2)= crsArray(2)*g(2);
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
                a(11,agn1)= crsArray(agn1)*g(agn1);
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
    % Change the value of cap to favor attachment of particular aggregate length  
    a(5,rst(5):re(5))= cap*g0*y(rst(5):re(5));
%     AggChangeLength = 10;
%     a(5,rst(5):AggChangeLength)= cap*g0*y(rst(5):AggChangeLength);
%     a(5,AggChangeLength+1:re(5))= cap*g0*y(AggChangeLength+1:re(5)).*(AggChangeLength+1:re(5));
%     a(5,AggChangeLength+1:re(5))= cap*g0*y(AggChangeLength+1:re(5))./(AggChangeLength+1:re(5));
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


F = F*vScale/v;
M = M*vScale/v;
y0t = y0t*vScale/v;
g0t = g0t*vScale/v;

% Real world time: treal = t(m)*v/3600;
t = t/3600/24;      % Covert Gillespie time to real time 
treal = t(end);
xFactor = floor(log10(t(end)));
t = t/(10^xFactor);