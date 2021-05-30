function [ Tp ] = assemble_temperature( ENL , NL )

% We want to assemble the prescribed  temp, similar to how we
% prescribe position in Mechanical problem

PD = size(NL,2);
NoN = size(NL,1);
DOCs = 0;

for i=1:NoN
    
    if (ENL (i,3) == -1) %if we have a neuman BC
        DOCs = DOCs + 1; %count DOFs
        Tp(DOCs,1) = ENL(i,PD+4);  %assign to this output vector  the Dirichlet value from ENL
        
    end
    
    
end




end