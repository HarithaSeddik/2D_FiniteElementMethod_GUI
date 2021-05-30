function [ Qp ] = assemble_flux( ENL , NL )
% We want to assemble the prescribed  heat fluxes, similar to how we
% prescribe forces in Mechanical problem

PD = size(NL,2);
NoN = size(NL,1);
DOFs = 0;

for i=1:NoN
    
    if (ENL (i,3) == 1) %if we have a neuman BC
        DOFs = DOFs + 1; %count DOFs
        
        Qp(DOFs,1) = ENL(i,PD+5);  %assign to this output vector  the neuman value from ENL
        
    end
    
    
end



end