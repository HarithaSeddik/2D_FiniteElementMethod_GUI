%This is an appdesigner function that determines which meshing function we
%call

function [NL, EL] = mesh_it(a1, b1, c1, d1, e1, elementType)

if(elementType == "2DQU4N")
    [NL, EL] = mesh_mk(a1, b1, c1, d1, e1);
elseif(elementType == "2DTR3N")
    [NL, EL] = mesh_3(a1, b1, c1, d1, e1);
elseif(elementType == "2DTR6N")
    [NL, EL] = mesh_6(a1, b1, c1, d1, e1);
elseif(elementType == "2DQU8N")
    [NL, EL] = mesh_8(a1, b1, c1, d1, e1);
elseif(elementType == "2DQU9N")
    [NL, EL] = mesh_mk_9(a1, b1, c1, d1, e1);
end

end
