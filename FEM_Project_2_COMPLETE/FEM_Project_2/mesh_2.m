function  [ NL, EL ] = mesh_2 (d1,d2,p,m ,R)

q = [ 0.5 - d1/2 0.5 - d2/2 ; 0.5 + d1/2 0.5 - d2/2 ; ...
      0.5 - d1/2 0.5 + d2/2 ; 0.5 + d1/2 0.5 + d2/2]; %corners of our domain

PD = 2;

NoN = 2*(p+1)*(m+1) + 2*(p-1)*(m+1);  % 2 partitions(4) since 2 of them will count the edges and the other 2 will not 

NoE = 4 * p * m;

NPE = 4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Nodes    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NL = zeros(NoN,PD);

 a = (q(2,1) - q(1,1)) / p;
 b = (q(3,2) - q(1,2)) / p;
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % NL of region 1 (bottom)
 
 coor11 = zeros ( (p+1)*(m+1) , PD);  
 
 for i = 1: p+1 %% nodes at the bottom of region 1
     
     coor11(i,1) = q(1,1) + (i-1) * a;
     coor11(i,2) = q(1,2) ;
 end
 
 for i= 1: p+1   %%% nodes on the circle end
     
     coor11(m*(p+1)+i,1) = R*cos( (5*pi/4) + (i-1)*(pi/2)/p) + 0.5 ; 
     coor11(m*(p+1)+i,2) = R*sin( (5*pi/4) + (i-1)*(pi/2)/p) + 0.5 ;
 
 end
 % Everything else in between
 for i= 1 : m-1
     
     for j = 1:p+1
         
         x = ( coor11(m*(p+1) + j , 1) - coor11(j,1) ) / m; %increment on along the edge in x
         y = ( coor11(m*(p+1) + j , 2) - coor11(j,2) ) / m;
         
         coor11(i*(p+1)+j ,1) = coor11((i-1)*(p+1)+j,1) + x;
         coor11(i*(p+1)+j ,2) = coor11((i-1)*(p+1)+j,2) + y;
     
     end
 end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % NL of region 2 (top)
 
 coor22 = zeros ( (p+1)*(m+1) , PD);  
 
 for i = 1: p+1 %% nodes at the bottom of region 1
     
     coor22(i,1) = q(3,1) + (i-1) * a;
     coor22(i,2) = q(3,2);
 end
 
 for i= 1: p+1   %%% nodes on the circle end
     
     coor22(m*(p+1)+i,1) = R*cos( (3*pi/4) - (i-1)*(pi/2)/p) + 0.5 ; 
     coor22(m*(p+1)+i,2) = R*sin( (3*pi/4) - (i-1)*(pi/2)/p) + 0.5 ;
 
 end
 
 for i= 1 : m-1
     
     for j = 1:p+1
         
         x = ( coor22(m*(p+1) + j , 1) - coor22(j,1) ) / m; %increment on along the edge in x
         y = ( coor22(m*(p+1) + j , 2) - coor22(j,2) ) / m;
         
         coor22(i*(p+1)+j ,1) = coor22((i-1)*(p+1)+j,1) + x;
         coor22(i*(p+1)+j ,2) = coor22((i-1)*(p+1)+j,2) + y;
     
     end
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % NL of region 3 (left)
 
 coor33 = zeros ( (p-1)*(m+1) , PD);  
 
 for i = 1: p-1 %% nodes at the bottom of region 3
     
     coor33(i,1) = q(1,1);
     coor33(i,2) = q(1,2) + (i)*b;
 end
 
 
 for i= 1: p-1   %%% nodes on the circle end
     
     coor33(m*(p-1)+i,1) = R*cos( (5*pi/4) - (i)*(pi/2)/p) + 0.5 ; 
     coor33(m*(p-1)+i,2) = R*sin( (5*pi/4) - (i)*(pi/2)/p) + 0.5 ;
 
 end
 
 for i= 1 : m-1
     
     for j = 1:p-1
         
         x = ( coor33(m*(p-1) + j , 1) - coor33(j,1) ) / m; %increment on along the edge in x
         y = ( coor33(m*(p-1) + j , 2) - coor33(j,2) ) / m;
         
         coor33(i*(p-1)+j ,1) = coor33((i-1)*(p-1)+j,1) + x;
         coor33(i*(p-1)+j ,2) = coor33((i-1)*(p-1)+j,2) + y;
     
     end
 end
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % NL of region 4 (right)
 
 coor44 = zeros ( (p-1)*(m+1) , PD);  
 
 for i = 1: p-1 %% nodes at the bottom of region 4
     
     coor44(i,1) = q(2,1);
     coor44(i,2) = q(2,2) + (i)*b;
 end
 
 
 for i= 1: p-1   %%% nodes on the circle end
     
     coor44(m*(p-1)+i,1) = R*cos( (7*pi/4) + (i)*(pi/2)/p) + 0.5 ; 
     coor44(m*(p-1)+i,2) = R*sin( (7*pi/4) + (i)*(pi/2)/p) + 0.5 ;
 
 end
 
 for i= 1 : m-1
     
     for j = 1:p-1
         
         x = ( coor44(m*(p-1) + j , 1) - coor44(j,1) ) / m; %increment on along the edge in x
         y = ( coor44(m*(p-1) + j , 2) - coor44(j,2) ) / m;
         
         coor44(i*(p-1)+j ,1) = coor44((i-1)*(p-1)+j,1) + x;
         coor44(i*(p-1)+j ,2) = coor44((i-1)*(p-1)+j,2) + y;
     
     end
 end
 
 %Now we store all these matrices together in one node list. 
 
 %NL = [ coor11; coor22; coor33; coor44];
 
 %This is okay but very difficult to generate a logical node list. 
 
 %This for loop arranges our node list according to the 'strips' method
 for i=1:m+1
     
     NL((i-1)*4*p+1:(i)*4*p,:) = [coor11((i-1)*(p+1) + 1      : (i)*(p+1),:);
                                  coor44((i-1)*(p-1) + 1      : (i)*(p-1),:);
                                  flipud(coor22((i-1)*(p+1)+1 : (i)*(p+1),:));
                                  flipud(coor33((i-1)*(p-1)+1 : (i)*(p-1),:))];   %flipud: flip up to down (first row to last row, last row to first row)
                              
 end
 
                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Elements    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EL = zeros(NoE,NPE);

%Strategy for Naming:
%(1) We go over each strip(belt) in the matrix
%(2) In each block we have local node/element numbers where we start counting 
% counter clockwise 
%(3) then we assign global node numbers from counter clockwise of each
%strip 
    %That is, if p = 4, and we have 4 regions, then the last node of the
    %strip has the global number 4*p (16)
%(4) For the next strip, we strt with the

for i= 1:m
    
    for j= 1: 4*p
        
        if j==1  %the first element in strip
            
            EL((i-1)*(4*p)+j,1) = (i-1)*(4*p)+j;   %first node is a multiple of 4*p, but excluding the first one
            EL((i-1)*(4*p)+j,2) = EL((i-1)*(4*p)+j,1) + 1    ; %2nd node is first node + 1
            EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j,1) + 4*(p)   ; %4th node is 1st node + 4*p
            EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,4) + 1    ; %3rd node is 4th node + 1
        
        elseif j==4*p %if last element in strip
            %note that we use i not i-1 since we want the contribution of 4*p from the get go
            EL((i-1)*(4*p)+j,1) = (i-1)*(4*p) + 1;   %first element is last row's first element we could say
            EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j,1) + 4*p - 1 ; %4th element is 1'st element + 1 row - 1 element
            EL((i-1)*(4*p)+j,2) = EL((i-1)*(4*p)+j,4) + 1 ; %2nd element is 4th element + 1
            EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,4)+ 4*p ;%  3rd element is 4th element + 1 row (4*p)
            
             
        else
            
            EL((i-1)*(4*p)+j,1) = EL((i-1)*(4*p)+j-1,2)  ;  %first node of current element is 2nd node from prev element
            EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j-1,3)  ;  % 4th node of current element is 3rd element from prev element
            EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,4) + 1;  %3rd node of current element is 4th node of current element + 1
            EL((i-1)*(4*p)+j,2) = EL((i-1)*(4*p)+j,1) + 1;  %2nd node of curren element is 1st node of current element + 1
            
        end
            
    end
end


 end

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 