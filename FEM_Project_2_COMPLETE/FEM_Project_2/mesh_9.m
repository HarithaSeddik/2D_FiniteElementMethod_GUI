%This is for our circle

function [NL, EL] = mesh_9(d1,d2,p,m, R)

%previously q was the corners of our domain
%now we use q to fix the circle (which is in the middle)
q = [0.5-d1/2 0.5-d2/2; 0.5+d1/2 0.5-d2/2; 0.5-d1/2 0.5+d2/2; 0.5+d1/2 0.5+d2/2];

PD = 2;

NoN = 2*(2*p+1)*(2*m+1) + 2*(2*p-1)*(2*m+1);
NoE = 4*p*m;
NPE = 9;

% Nodes
NL = zeros(NoN, PD);

a = (q(2,1) - q(1,1))/(2*p);
b = (q(3,2) - q(1,2))/(2*p);

%Region 1
coor11 = zeros((2*p+1)*(2*m+1), PD);

%The flat part of region 1 at the bottom
for i = 1:2*p+1
    coor11(i, 1) = q(1,1) + (i-1)*a;
    coor11(i, 2) = q(1,2);
end

%the innermost curved part of region 1 at the top
for i = 1:2*p+1
    coor11(2*m*(2*p+1) + i, 1) = R*cos( (5*pi/4) + (i-1)*(pi/2)/(2*p) ) + 0.5;
    coor11(2*m*(2*p+1) + i, 2) = R*sin( (5*pi/4) + (i-1)*(pi/2)/(2*p) ) + 0.5;
end

%everything else in between
for i = 1:2*m-1
    for j = 1:2*p+1
        
        x = (coor11(2*m*(2*p+1) + j, 1) - coor11(j, 1))/(2*m);
        y = (coor11(2*m*(2*p+1) + j, 2) - coor11(j, 2))/(2*m);
        
        coor11(i*(2*p+1)+j, 1) = coor11((i-1)*(2*p+1)+j, 1) + x;
        coor11(i*(2*p+1)+j, 2) = coor11((i-1)*(2*p+1)+j, 2) + y;
    end
end


%Region 2
coor22 = zeros((2*p+1)*(2*m+1), PD);

%The flat part of region 1 at the bottom
for i = 1:2*p+1
    coor22(i, 1) = q(3,1) + (i-1)*a;
    coor22(i, 2) = q(3,2);
end

%the innermost curved part of region 1 at the top
for i = 1:2*p+1
    coor22(2*m*(2*p+1) + i, 1) = R*cos( (3*pi/4) - (i-1)*(pi/2)/(2*p) ) + 0.5;
    coor22(2*m*(2*p+1) + i, 2) = R*sin( (3*pi/4) - (i-1)*(pi/2)/(2*p) ) + 0.5;
end

%everything else in between
for i = 1:2*m-1
    for j = 1:2*p+1
        
        x = (coor22(2*m*(2*p+1) + j, 1) - coor22(j, 1))/(2*m);
        y = (coor22(2*m*(2*p+1) + j, 2) - coor22(j, 2))/(2*m);
        
        coor22(i*(2*p+1)+j, 1) = coor22((i-1)*(2*p+1)+j, 1) + x;
        coor22(i*(2*p+1)+j, 2) = coor22((i-1)*(2*p+1)+j, 2) + y;
    end
end

%Region 3
coor33 = zeros((2*p+1)*(2*m+1), PD);

%The flat part of region 1 at the bottom
for i = 1:2*p-1
    coor33(i, 1) = q(1,1);
    coor33(i, 2) = q(1,2) + i*b;
end

%the innermost curved part of region 1 at the top
for i = 1:2*p-1
    coor33(2*m*(2*p-1) + i, 1) = R*cos( (5*pi/4) - (i)*(pi/2)/(2*p) ) + 0.5;
    coor33(2*m*(2*p-1) + i, 2) = R*sin( (5*pi/4) - (i)*(pi/2)/(2*p) ) + 0.5;
end

%everything else in between
for i = 1:2*m-1
    for j = 1:2*p-1
        
        x = (coor33(2*m*(2*p-1) + j, 1) - coor33(j, 1))/(2*m);
        y = (coor33(2*m*(2*p-1) + j, 2) - coor33(j, 2))/(2*m);
        
        coor33(i*(2*p-1)+j, 1) = coor33((i-1)*(2*p-1)+j, 1) + x;
        coor33(i*(2*p-1)+j, 2) = coor33((i-1)*(2*p-1)+j, 2) + y;
    end
end


%Region 4
coor44 = zeros((2*p-1)*(2*m+1), PD);

%The flat part of region 1 at the bottom
for i = 1:2*p-1
    coor44(i, 1) = q(2,1);
    coor44(i, 2) = q(2,2) + i*b;
end

%the innermost curved part of region 1 at the top
for i = 1:2*p-1
    coor44(2*m*(2*p-1) + i, 1) = R*cos( (7*pi/4) + (i)*(pi/2)/(2*p) ) + 0.5;
    coor44(2*m*(2*p-1) + i, 2) = R*sin( (7*pi/4) + (i)*(pi/2)/(2*p) ) + 0.5;
end

%everything else in between
for i = 1:2*m-1
    for j = 1:2*p-1
        
        x = (coor44(2*m*(2*p-1) + j, 1) - coor44(j, 1))/(2*m);
        y = (coor44(2*m*(2*p-1) + j, 2) - coor44(j, 2))/(2*m);
        
        coor44(i*(2*p-1)+j, 1) = coor44((i-1)*(2*p-1)+j, 1) + x;
        coor44(i*(2*p-1)+j, 2) = coor44((i-1)*(2*p-1)+j, 2) + y;
    end
end



%%%Now we store everything in a node list

%This is okay, but complicated to generate a element list
%NL = [coor11; coor22; coor33; coor44];

%We reorder it to make generating element list easier
for i = 1:2*m+1
    NL((i-1)*8*p+1:(i)*8*p, :) = [coor11((i-1)*(2*p+1)+1:(i)*(2*p+1),:);
                                  coor44((i-1)*(2*p-1)+1:(i)*(2*p-1),:);
                                  flipud(coor22((i-1)*(2*p+1)+1:(i)*(2*p+1),:));
                                  flipud(coor33((i-1)*(2*p-1)+1:(i)*(2*p-1),:))];
end



% Elements

EL = zeros(NoE,NPE);

%The strategy for naming here is:
%1. We go over each strip in the matrix
%2. Each block we have local node numbers that start and count
%counterclockwise (1, 2, 3, 4)
%3. We then assign global node numbers from counterclockwise of each strip
    %That is, if p = 4, and we have 4 regions, then the last node of the
    %strip has a global number of 16
%4. For the next strip, we start w ith the next number (Here it's 17)
    %That means that for the first block:
    %Local Node 1 = Global Node 1
    %Local Node 2 = Global Node 12
    %Local Node 3 = Global Node 17
    %Local Node 4 = Global Node 18
    


for i= 1:m
    
    for j= 1: 4*p
        
        if j==1  %the first element in strip
            
            EL((i-1)*(4*p)+j,1) = (i-1)*(16*p)+1;
            EL((i-1)*(4*p)+j,2) = (i-1)*(16*p)+3;
            EL((i-1)*(4*p)+j,3) = (i)*(16*p)+3;
            EL((i-1)*(4*p)+j,4) = (i)*(16*p)+1;
            EL((i-1)*(4*p)+j,5) = (i-1)*(16*p)+2;
            EL((i-1)*(4*p)+j,6) = (i)*(8*p)+3;
            EL((i-1)*(4*p)+j,7) = (i)*(16*p)+2;
            EL((i-1)*(4*p)+j,8) = (i)*(8*p)+1;
            EL((i-1)*(4*p)+j,9) = (i-1)*(16*p) + 8*p + 2;
        
        elseif j==4*p %if last element in strip
            
            EL((i-1)*(4*p)+j,1) = (i-1)*(16*p) + (8*p) - 1 ;   
            EL((i-1)*(4*p)+j,2) = EL((i-1)*(4*p)+j,1) - (8*p) +2 ; 
            EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,2) + 16*p;
            EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j,1) + 16*p;
            
            EL((i-1)*(4*p)+j,5) = EL((i-1)*(4*p)+j,1) + 1;
            EL((i-1)*(4*p)+j,6) = EL((i-1)*(4*p)+j,2) + 8*p;
            EL((i-1)*(4*p)+j,7) = EL((i-1)*(4*p)+j,4) + 1;
            EL((i-1)*(4*p)+j,8) = EL((i-1)*(4*p)+j,1) + 8*p;
            EL((i-1)*(4*p)+j,9) = EL((i-1)*(4*p)+j,1) + (8*p)+1;
            
             
        else
            
            EL((i-1)*(4*p)+j,1) = EL((i-1)*(4*p)+j-1,1) +2;  
            EL((i-1)*(4*p)+j,2) = EL((i-1)*(4*p)+j-1,2) +2;  
            EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j-1,3) +2;  
            EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j-1,4) +2;
            EL((i-1)*(4*p)+j,5) = EL((i-1)*(4*p)+j-1,5) +2;
            EL((i-1)*(4*p)+j,6) = EL((i-1)*(4*p)+j-1,6) +2;
            EL((i-1)*(4*p)+j,7) = EL((i-1)*(4*p)+j-1,7) +2;
            EL((i-1)*(4*p)+j,8) = EL((i-1)*(4*p)+j-1,8) +2;
            EL((i-1)*(4*p)+j,9) = EL((i-1)*(4*p)+j-1,9) +2;
            
        end
            
    end
end


end

