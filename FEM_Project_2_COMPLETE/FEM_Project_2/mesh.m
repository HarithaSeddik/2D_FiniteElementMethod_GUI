function  [ NL, EL ] = mesh( d1,d2,p,m )

q = [0 0; d1 0; 0 d2 ; d1 d2];

NoN = (p+1) * (m+1); %Number of Nodes

NoE = p * m; %Number of Elements

NPE = 2; %Nodes per element, 4 since we have a quadratic 

PD = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Nodes    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NL = zeros(NoN,NPE);  

a = (q(2,1) - q(1,1))/p; % horiz length between nodes
b = (q(3,2) - q(1,2))/m; % vertcal length between nodes

for i = 1: m+1  %counting the rows
    
    for j = 1: p+1 %counting the colums
        
        NL((i-1) * (p+1) + j,1) = q(1,1) + (j-1)*a;
        NL((i-1) * (p+1) + j,2) = q(1,2) + (i-1)*b;
        
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Elements    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EL = zeros(NoE, NPE );

for i = 1:m
    for j = 1:p
        
        if j == 1 %if we are at the beginning of each column, identify the nodes
                    
                EL((i-1)*p + j,1) = (i-1)*(p+1) + j; %bottom left corner
                EL((i-1)*p + j,2) = EL((i-1)*p + j,1) + 1;  %bottom right corner
                EL((i-1)*p + j,4) = EL((i-1)*p + j,1) + p+1;  %top left corner
                EL((i-1)*p + j,3) = EL((i-1)*p + j,4) + 1;  %top right corner
                
        else 
                EL((i-1)*p + j,1) = EL((i-1)*p+j -1,2); %first node of this element is 2nd node of previous element
                EL((i-1)*p + j,4) = EL((i-1)*p+j -1,3); %fourth node of this element is the 3rd node of previous elemet
                EL((i-1)*p + j,2) = EL((i-1)*p + j,1) + 1;
                EL((i-1)*p + j,3) = EL((i-1)*p + j,4) + 1;
       
        end
    end



end
end
