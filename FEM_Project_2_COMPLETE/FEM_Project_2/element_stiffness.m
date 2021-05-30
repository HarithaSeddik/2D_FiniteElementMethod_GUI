function output = element_stiffness( nl,NL)  

%Function that calculates the stiffness of single 2D finite element

% nl contains the nodes connected to that 2D element

    NPE = size(nl,2);  %nodes per element
    PD = size(NL,2);  %problem dimension
    
    %coor = x';  %for javilis convention
    
    V = getGlobalv;
    GPE = getGlobalgpe;
    E = getGlobalE;
    

    
    coor = zeros(PD,NPE); %coordinates of each element
    for i = 1:NPE
        for j = 1:PD
            coor(j,i) = NL(nl(i),j);
        end
    end
    
    
switch NPE
    case 3
        GPE = 1;
    case 4
        GPE = 4;
    case 6
        GPE = 3;
    case 8
        GPE = 4;
    case 9
        GPE = 9;
end


    
    k = zeros(PD*NPE,PD*NPE); %initialize element stiffness
    
    for i= 1:NPE %loop over nodes
        for j= 1:NPE  %loop over nodes
           S = zeros(PD, PD); %%or rather, K
            for gp = 1:GPE %loop over Gauss Points
                
                J = zeros(PD,PD);  %JACOBIAN
                grad = zeros(PD,NPE);
                
                [xi,eta,alpha] = GaussPoint(NPE,GPE,gp);
                grad_nat = grad_N_nat( NPE, xi , eta);
                J = coor * grad_nat';

                grad = inv(J)' * grad_nat;
               
               for a = 1:PD
                 for c = 1:PD
                       for b = 1:PD
                          for d = 1:PD
                             S(a,c) = S(a,c) + grad(b, i) * constitutive(a, b, c, d,E,V) * grad(d, j) * det(J) * alpha;
                          end
                       end
                 end
               end
            
            end
 
            k((i-1)*PD+1:i*PD, (j-1)*PD+1:j*PD) = S;
         end
    output = k;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = grad_N_nat( NPE, xi , eta)
%returns matrix representing the gradient of the shape functions in natural
%domain evaluated at xi and eta

    PD = 2; 
    
    result = zeros(PD,NPE);
    
    if (NPE ==3) %triangular
        
        result(1,1) = 1;
        result(1,2) = 0;
        result(1,3) = -1;
        
        result(2,1) = 0;
        result(2,2) = 1;
        result(2,3) = -1;
        
    elseif(NPE == 4)
        
        result(1,1) = -1/4*(1-eta);
        result(1,2) =  1/4*(1-eta);
        result(1,3) =  1/4*(1+eta);
        result(1,4) = -1/4*(1+eta);
        
        result(2,1) = -1/4*(1-xi);
        result(2,2) = -1/4*(1+xi);
        result(2,3) =  1/4*(1+xi);
        result(2,4) =  1/4*(1-xi);
        
    elseif(NPE == 6)
        
       result(1,1)=-1 +(4*xi);
       result(1,2)=0;
       result(1,3)=-3+(4*xi)+(4*eta);
       result(1,4)=4*eta;
       result(1,5)=-4*eta;
       result (1,6)= -4*(-1+eta+(2*xi));
       
       result(2,1)=0;
       result(2,2)=-1 +(4*eta);
       result(2,3)=-3+(4*xi)+(4*eta);
       result(2,4)=4*xi;
       result(2,5)=-4*(-1+xi+(2*eta));
       result(2,6)= -4*xi;
        
    elseif (NPE == 9)
        
        result(1,1) = (1/4)*(1-2*xi)*(1-eta)*eta;
        result(1,2) = -(1/4)*(1+2*xi)*(1-eta)*eta;
        result(1,3) = (1/4)*(1+2*xi)*(1+eta)*eta;
        result(1,4) = -(1/4)*(1-2*xi)*(1+eta)*eta;
        result(1,5) = xi*(1-eta)*eta;
        result(1,6) = (1/2)*(1+2*xi)*(1-eta)*(1+eta);
        result(1,7) = -xi*(1-eta)*eta;
        result(1,8) = -(1/2)*(1-2*xi)*(1-eta)*(1+eta);
        result(1,9) = -2*xi*(1-eta)*(1+eta);
        
        result(2,1) = (1/4)*(1-2*eta)*(1-xi)*xi;
        result(2,2) = -(1/4)*(1-2*eta)*(1+xi)*xi;
        result(2,3) = (1/4)*(1+2*eta)*(1+xi)*xi;
        result(2,4) = -(1/4)*(1+2*eta)*(1-xi)*xi;
        result(2,6) = -xi*(1+xi)*eta;
        result(2,5) = (1/2)*(-1+2*eta)*(1-xi)*(1+xi);
        result(2,8) = -eta*(1-eta)*eta;
        result(2,7) = (1/2)*(1+2*eta)*(1-xi)*(1+xi);
        result(2,9) = -2*eta*(1-xi)*(1+xi);
        
        
    elseif (NPE == 8)
        
        result(1,1) = (1/4)*(2*xi+eta)*(1-eta);
        result(1,2) = (1/4)*(2*xi-eta)*(1-eta);
        result(1,3) = (1/4)*(2*xi+eta)*(1+eta);
        result(1,4) = (1/4)*(2*xi-eta)*(1+eta);
        result(1,5) = -xi*(1-eta);
        result(1,6) = (1/2)*(1+eta)*(1-eta);
        result(1,7) = -xi*(1+eta);
        result(1,8) = (-1/2)*(1+eta)*(1-eta);
        
        result(2,1) = (1/4)*(2*eta+xi)*(1-xi);
        result(2,2) = (1/4)*(2*eta-xi)*(1+xi);
        result(2,3) = (1/4)*(2*eta+xi)*(1+xi);
        result(2,4) = (1/4)*(2*eta-xi)*(1-xi);
        result(2,6) = -eta*(1+xi);
        result(2,5) = -(1/2)*(1+xi)*(1-xi);
        result(2,8) = -eta*(1-xi);
        result(2,7) = (1/2)*(1+xi)*(1-xi);

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ xi , eta, alpha ] = GaussPoint( NPE , GPE , gp)

        if (NPE == 3 || NPE == 6)  %D2TR3N
            
            if GPE == 1
            
            xi = 1/3; eta = 1/3 ; alpha = 1 * (1/2) ; 
            
            elseif GPE == 3
                
                if gp == 1
                    xi = 1/6; eta = 1/6 ; alpha = (1/3) * (1/2) ;
                elseif gp == 2
                    xi = 2/3; eta = 1/6 ; alpha = (1/3) * (1/2) ;
                elseif gp == 3
                    xi = 1/6; eta = 2/3 ; alpha = (1/3) * (1/2) ;
                end
            end
        
        elseif (NPE == 4 || NPE == 8)  %D2QU4N (or any D2QU
            
            if (GPE == 1) %if 1 gauss points
                
                xi = 0 ; eta = 0; alpha = 4;
                
            elseif (GPE == 4) %if 4 gauss points
                
                if (gp == 1)  %first gp
                    
                    xi = -1/sqrt(3); eta= -1/sqrt(3); alpha = 1;
                    
                elseif (gp == 2)
                    
                    xi =  1/sqrt(3); eta= -1/sqrt(3); alpha = 1;
                    
                elseif (gp == 3)
                    
                    xi =  1/sqrt(3); eta= 1/sqrt(3); alpha = 1;
                elseif (gp == 4)
                    
                    xi = -1/sqrt(3); eta= 1/sqrt(3); alpha = 1;
                end
            end
            
                  
        elseif NPE==9
            if (GPE == 1) %if 1 gauss points
                xi = 0 ; eta = 0; alpha = 4;
            elseif (GPE == 9)
                if gp==1
                    xi = -sqrt(3/5); eta = -sqrt(3/5); alpha = (5/9)*(5/9);
                elseif gp==2
                    xi = +sqrt(3/5); eta = -sqrt(3/5); alpha = (5/9)*(5/9);
                elseif gp==3
                    xi = +sqrt(3/5); eta = sqrt(3/5); alpha = (5/9)*(5/9);
                elseif gp==4
                    xi = -sqrt(3/5); eta = +sqrt(3/5); alpha = (5/9)*(5/9);
                elseif gp==5
                    xi = 0; eta = -sqrt(3/5); alpha = (5/9)*(8/9);
                elseif gp==6
                    xi = +sqrt(3/5); eta =0; alpha = (5/9)*(8/9);
                elseif gp==7
                    xi = 0; eta = sqrt(3/5); alpha = (5/9)*(8/9);
                elseif gp==8
                    xi = -sqrt(3/5); eta = 0; alpha = (5/9)*(8/9);
                elseif gp==9
                    xi = 0; eta = 0; alpha = (5/9)*(8/9);    
                end
            end
        end
end

function C = constitutive(i, j, k, l ,E, V)
% E = 80/3;
% nu = 1/3;
C = (E/(2*(1+V))) * ((i==l)*(j==k) + (i==k)*(j==l)) + E*V/(1-V^2) * (i==j)*(k==l);
end

function v = getGlobalv
global V
v = V;
end

function GPE = getGlobalgpe
global gpe
GPE = gpe;
end

function e = getGlobalE
global E
e = E;
end