function [NL,EL] = mesh_8(d1,d2,p,m,R) 
tic
COC = [d1/2 d2/2] ;
m = 2*m ;
%assembling the NL 
o =1 ;
REF = [0 0];

%defining the polygon in the center:
POLY = nsidedpoly(4*p, 'Center', COC, 'SideLength', 2*R*sin(pi/(4*p)) ) ;

V = POLY.Vertices; 

VROT = zeros(size(V));
 
if floor(p/2)==p/2 
    % rotate the vertices if p is even 
for i = 1:size(V,1)
       
    VROT(i,1) = (V(i,1)-COC(1,1)) * cos(pi/(4*p)) - (V(i,2) - COC(1,2)) *sin(pi/(4*p)) + COC(1,1);
    VROT(i,2) = (V(i,1)-COC(1,1)) * sin(pi/(4*p)) + (V(i,2) - COC(1,2)) *cos(pi/(4*p)) + COC(1,2);
    
end 
V = VROT ;
end
%polygon in the center is ready, now lets rearrange the coordinates in the
%matrix V to fit our model

if floor(p/2)==p/2 
    a = p+1;
    S = ceil(a/2) ;
else
    S = ceil(p/2) ;
end    

if p == 1
   V(2:end,:) = flipud(V(2:end,:)) ;
else 
    V(1:S,:) = flipud(V(1:S,:)) ;
    V(S+1:end, :) = flipud(V(S+1:end, :)) ;
end


%from 0 to the side
%from 0 to the side
for i = 1 : p 
    
    NL(i,:) = REF ;
    
    for j = 1:m
        x = REF(1) + (V(i,1) - REF(1))*(j/m) ;
        y = REF(2) + (V(i,2) - REF(2))*(j/m) ;
        NL(i+4*p*j, :) = [x y] ;
    end
    
    REF(1,1) = REF(1,1) + d1/p ;
end
%from the lower right corner upwards
for i = 1+p : 2*p 
    
    NL(i,:) = REF ;
    
    for j = 1:m
        x = REF(1) + (V(i,1) - REF(1))*(j/m) ;
        y = REF(2) + (V(i,2) - REF(2))*(j/m) ;
        NL(i+4*p*j, :) = [x y] ;
    end
    
    REF(1,2) = REF(1,2) + d2/p ;
end
%from right upper corner to the left corner
for i = 2*p+1 : 3*p 
    
    NL(i,:) = REF ;
    
    for j = 1:m
        x = REF(1) + (V(i,1) - REF(1))*(j/m) ;
        y = REF(2) + (V(i,2) - REF(2))*(j/m) ;
        NL(i+4*p*j, :) = [x y] ;
    end
    
    REF(1,1) = REF(1,1) - d1/p ;
end
%from the left upper corner down
for i = 3*p+1 : 4*p 
    
    NL(i,:) = REF ;
    
    for j = 1:m
        x = REF(1) + (V(i,1) - REF(1))*(j/m) ;
        y = REF(2) + (V(i,2) - REF(2))*(j/m) ;
        NL(i+4*p*j, :) = [x y] ;
    end
    
    REF(1,2) = REF(1,2) - d2/p ;
end

OP = [] ;
c = 1;
for i = 1 : 2*4*p : (4*p*m + 4*p)
   for j = i : (c*4*p-1)
    x = NL(j,1) + (NL(j+1,1) - NL(j,1))*(0.5) ;
    y = NL(j,2) + (NL(j+1,2) - NL(j,2))*(0.5) ;
    OP = [OP ; x y];
    if j == (c*4*p-1)
        x = NL(j+1,1) + (NL(j+1-(4*p-1),1) - NL(j+1,1))*(0.5) ;
        y = NL(j+1,2) + (NL(j+1-(4*p-1),2) - NL(j+1,2))*(0.5) ;
        OP = [OP ; x y];
    end
   end
    c = c+2; 
end
NL = [NL; OP] ;

%make EL
 m = m /2 ;
EL = [o o+1 (2*4*p+1+o) (2*4*p+o) (2*4*p*m + 4*p + o) (o+4*p+1) (2*4*p*m + 2*4*p + o) (o+4*p)] ;
j = 1;
for i = 2 : 4*p*m
   EL(i,:) = EL(i-1,:) + o ;%[EL(1)+o EL(2)+o EL(3)+o EL(4)+o EL(5)+o EL(6)+o EL(7)+o EL(8)+o];
   
   for k = 1: size(EL,1)
        if EL(k,1) == 4*p*j
       EL(k, 2) = EL(k-1,2) - (4*p-1) ; 
       EL(k, 3) = EL(k-1,3) - (4*p-1) ; 
       EL(k, 6) = EL(k-1,6) - (4*p-1) ; 
        end
   if EL(k,1) == 4*p*j+1
       
       EL(k,1) = EL(k-1,1) + 4*p+1 ;
       EL(k,2) = EL(k-1,2) + 2*4*p+1 ;
       EL(k,3) = EL(k-1,3) + 2*4*p+1 ;
       EL(k,4) = EL(k-1,4) + 4*p+1 ;
       EL(k,6) = EL(k-1,6) + 2*4*p+1 ;
       EL(k,8) = EL(k-1,8) + 4*p+1 ;
       
       j = j + 2;
   end
   end
   
end
toc
end