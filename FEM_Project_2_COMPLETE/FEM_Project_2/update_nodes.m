function [ ENL ] = update_nodes( ENL , U , NL )

PD = size(NL, 2);
NoN = size(NL, 1);
DOFs = 0;
DOCs = 0;

for i = 1:NoN
    for j = 1:PD
        if(ENL(i, j + PD) == 1)
            DOFs = DOFs + 1;
            ENL(i, j + PD*4) = U(DOFs);
        end
    end
end


end