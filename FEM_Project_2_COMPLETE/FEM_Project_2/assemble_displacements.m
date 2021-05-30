function [ Up ] = assemble_displacements( ENL , NL )

NoN = size(NL, 1);
PD = size(NL, 2);
DOCs = 0;

%only assembles displacements in x-direction

for i = 1:NoN
    for j = 1:PD
        if(ENL(i, j+PD) == -1)
            DOCs = DOCs + 1;
            Up(DOCs, 1) = ENL(i, j+4*PD);
        end
    end
end



end