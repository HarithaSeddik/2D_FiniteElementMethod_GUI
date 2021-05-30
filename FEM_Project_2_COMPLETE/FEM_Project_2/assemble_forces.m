function [ Fp ] = assemble_forces( ENL , NL )

NoN = size(NL, 1);
PD = size(NL, 2);
DOFs = 0;

%Only assembles forces in x-direction

for i = 1:NoN
    for j = 1:PD
        if(ENL(i,j+PD) == 1)
            DOFs = DOFs + 1;
            Fp(DOFs, 1) = ENL(i, j+5*PD);
        end
    end
end


end