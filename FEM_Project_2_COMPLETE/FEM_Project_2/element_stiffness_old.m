function output = element_stiffness( nl,NL)  

%Function that calculates the stiffness of single 2D finite element

% nl contains the nodes connected to that 2D element

    NoN = size(NL,1); %number of nodes
    %NoE = size(EL,1); %number of elements
    NPE = size(nl,2);  %nodes per element
    PD = size(NL,2);  %problem dimension
    
    conductivity = getGlobalk;
    GPE = getGlobalgpe;
    
    %coor = x';  %for javilis convention
    
    %conductivity = 100;
    
    
    
%     if ( NPE == 3 )
%         GPE = 1;  % Gauss per element
%     elseif ( NPE == 4)
%         GPE = 4;  % Gauss per element
%     end
    
    coor = zeros(PD,NPE); %coordinates of each element
    for i = 1:NPE
        for j = 1:PD
            coor(j,i) = NL(nl(i),j);
        end
    end
    
    k = zeros(NPE,NPE); %initialize element stiffness
    
    for i= 1:NPE %loop over nodes
       
        for j= 1:NPE  %loop over nodes
            
            K = 0; %reset stiffness to zero
            
            for gp = 1:GPE %loop over Gauss Points
                
                J = zeros(PD,PD);  %JACOBIAN
                
                grad = zeros(PD,NPE);
                
                [xi,eta,alpha] = GaussPoint(NPE,GPE,gp);
                
                grad_nat = grad_N_nat( NPE, xi , eta);
                
                J = coor * grad_nat';

                grad = J^(-1) * grad_nat;
                
%                 for n= 1:NPE
%                     grad(:,n) = grad(:,n) + inv(J)' * grad_nat(:,n);
%                 end
                
                for b= 1:PD   %loop over PD since we have xi and eta 
                   
                    K = K + (-conductivity) * grad(b,i) * grad(b,j) * det(J) * alpha;
               
                end
                
                
            end
            
            k(i,j) = K;
            
        end
        
    end
    
    output = k;

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
        
        

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ xi , eta, alpha ] = GaussPoint( NPE , GPE , gp)

        if (NPE == 3)  %D2TR3N
            
            xi = 1/3; eta = 1/3 ; alpha = 1 * (1/2) ; 
        
        elseif (NPE == 4)  %D2QU4N (or any D2QU
            
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


        end
        
end

function conductivity = getGlobalk
global k
conductivity = k;
end

function GPE = getGlobalgpe
global gpe
GPE = gpe;
end

