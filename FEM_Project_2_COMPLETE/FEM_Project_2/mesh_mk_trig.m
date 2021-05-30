%This is for our circle

function [NL, EL] = mesh_mk_trig(a,b,p,m, R)

NL = [1.5 1.5];
SEL = [0 0 0 0];
EL = [0 0 0];

node_array = [1.5 1.5 1];
node_array(1,:) = [];

d1 = 1;
d2 = 1;
theta_step = 90/p;
node_count = 0;

for theta = 45 : theta_step : 45+(360-theta_step)
    
    if theta >= 45 && theta <= 135
        r_final = (d2/2)/sind(theta);
    end
    if theta > 135 && theta < 225
        r_final = -(d1/2)/cosd(theta);
    end
    if theta >= 225 && theta <= 315
        r_final = -(d2/2)/sind(theta);
    end
    if theta > 315 && theta < 405
        r_final = (d1/2)/cosd(theta);
    end
    
        
    for r = R : ((r_final - R)/(m)) : r_final
        
        %node = [r*cosd(theta) r*sind(theta)];
        NL = [NL ; [r*cosd(theta) r*sind(theta)]];
        node_count = node_count + 1;
        node_a = [node_count [r*cosd(theta) r*sind(theta)]];
        node_array = [node_array ; node_a];
    
    end
        
end


NL(1,:) = [];
SEL(1,:) = [];

NL(:,1) = a*NL(:,1);
NL(:,2) = b*NL(:,2);

b = zeros(size(node_array,1) , 1);
node_array = [node_array b];

for i = 1:size(node_array,1)
    node_array(i,4) = rem(node_array(i,1) , m+1);   
end

num_of_elements = 4*p*m;
n_c = 0;

for e = 1 : num_of_elements - m
    
    n_c = n_c + 1;
    
    if node_array(n_c , 4) == 0
        n_c = n_c + 1;
    end
    el_1 = n_c    ;
    el_2 = n_c + 1;
    el_3 = n_c + (m+2);
    el_4 = n_c + (m+1);
    
    %element = [el_1 el_2 el_3 el_4];
    SEL = [SEL ; [el_1 el_2 el_3 el_4]];
    
end

n_c = n_c + 1;

for i = 1:m
    
    n_c = n_c + 1;
    
    el_1 = n_c    ;
    el_2 = n_c + 1;
    el_3 = i + 1;
    el_4 = i    ;
    
    SEL = [SEL ; [el_1 el_2 el_3 el_4]];
    
end

EL(1,:) = [];

for i = 1:size(SEL,1)
    EL = [EL ; [SEL(i,1) SEL(i,2) SEL(i,4)]];
    EL = [EL ; [SEL(i,2) SEL(i,3) SEL(i,4)]];

end