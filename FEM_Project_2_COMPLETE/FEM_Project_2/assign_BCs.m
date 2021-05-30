%Changes to make:
    %Make prescribed displacements an unknown you can change (Global var)


function [ ENL , DOFs , DOCs ] = assign_BCs( NL, BC_flag, def_mag )


NoN = size(NL, 1);
PD = size(NL, 2);

ENL = zeros(NoN, PD + 5);
ENL(:, 1:PD) = NL;

switch BC_flag
    case 'Extension'
        %We apply extension by fixing one side and pulling the other
        for i = 1:NoN
            if(ENL(i, 1) == 0) %x = 0 is fixed to 0 displacement
                ENL(i, 3) = -1; %Dirichlet
                ENL(i, 4) = -1; %Dirichlet
                ENL(i, 9) = 0; %Prescribed Displacement
                ENL(i, 10) = 0; %Prescribed Displacement
                
            elseif(ENL(i,1) == 1) %x = 1 is fixed to 0.1 displacement only in x direction
                ENL(i,3) = -1; %Dirichlet
                ENL(i, 4) = -1; %Dirichlet
                ENL(i, 9) = def_mag; %Prescribed Displacement
                ENL(i, 10) = 0; %Prescribed Displacement
                
            else
                ENL(i, 3) = 1; %Neumann
                ENL(i, 4) = 1; %Neumann
                ENL(i, 11) = 0; %Prescribed Force
                ENL(i, 12) = 0; %Prescribed Force
            end
        end
                
        
    case 'Shear'
            
        for i = 1:NoN
            if(ENL(i, 2) == 0) %x = 0 is fixed to 0 displacement
                ENL(i, 3) = -1; %Dirichlet
                ENL(i, 4) = -1; %Dirichlet
                ENL(i, 9) = 0; %Prescribed Displacement
                ENL(i, 10) = 0; %Prescribed Displacement
                
            elseif(ENL(i,2) == 1) %x = 1 is fixed to 0.1 displacement only in x direction
                ENL(i,3) = -1; %Dirichlet
                ENL(i, 4) = -1; %Dirichlet
                ENL(i, 9) = def_mag; %Prescribed Displacement
                ENL(i, 10) = 0; %Prescribed Displacement
                
            else
                ENL(i, 3) = 1; %Neumann
                ENL(i, 4) = 1; %Neumann
                ENL(i, 11) = 0; %Prescribed Force
                ENL(i, 12) = 0; %Prescribed Force
            end
        end
        
    case 'Expansion'
        for i = 1:NoN
            if(ENL(i, 1) == 0 || ENL(i, 1)==1 || ENL(i, 2) == 0 || ENL(i, 2)==1 ) %if boundaries
                ENL(i, 3) = -1; %Dirichlet
                ENL(i, 4) = -1; %Dirichlet
                ENL(i, 9) =  ENL(i,1) * def_mag; %Prescribed X Displacement
                ENL(i, 10) = ENL(i,2) * def_mag; %Prescribed Displacement
                
            else
                ENL(i, 3) = 1; %Neumann
                ENL(i, 4) = 1; %Neumann
                ENL(i, 11) = 0; %Prescribed Force
                ENL(i, 12) = 0; %Prescribed Force
            end
        end
end

DOCs = 0;
DOFs = 0;

%Intialize the column of temporary degrees
for i = 1:NoN
    for j = 1:PD
        if(ENL(i, j + PD) == 1) %tells us that its a degree of freedom
            DOFs = DOFs + 1;
            ENL(i, j + PD*2) = DOFs;
        else
            DOCs = DOCs - 1;
            ENL(i, j + PD*2) = DOCs;
        end
    end
end

%Initialize column of global degrees
DOCs = abs(DOCs);

for i = 1:NoN
    for j = 1:PD
        if(ENL(i, j + PD) < 0) %dirichlet condition
            ENL(i, j + PD * 3) = abs(ENL(i, j+ 2*PD)) + DOFs;
        else
            ENL(i, j + PD*3) = ENL(i, j+ 2*PD);
        end
    end
end
    


end