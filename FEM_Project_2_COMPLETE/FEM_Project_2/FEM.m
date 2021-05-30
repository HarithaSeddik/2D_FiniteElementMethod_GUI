function [NL, EL, ENL]= FEM (NL, EL, BC_flag ,E , V, GPE,def_mag)
%%%%%  Function that takes care of the Process of 2D FEM, regardless of the
%%%%%  element type. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


setGlobalV(V)
setGlobalgpe(GPE)
setGlobalE(E)

   
[ENL,DOFs,DOCs] = assign_BCs(NL,BC_flag,def_mag) ; 
ENL


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   PROCESS.   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = assemble_stiffness(ENL,EL,NL);

Fp = assemble_forces(ENL,NL);

Up = assemble_displacements(ENL,NL);

K_reduced = K(1:DOFs,1:DOFs);

F = Fp - K(1:DOFs,DOFs+1:(DOFs+DOCs)) * Up;

U = K_reduced \ F;

ENL = update_nodes(ENL,U,NL);
ENL;

Fu = K(DOFs+1:DOFs+DOCs,1:DOFs) * U + K(DOFs+1:(DOFs+DOCs),DOFs+1:(DOFs+DOCs)) * Up;


end

function setGlobalV(val)
    global V
    V = val;
    
end

function setGlobalgpe(val)
    global gpe
    gpe = val;
    
end

function setGlobalE(val)
    global E
    E = val;
    
end
