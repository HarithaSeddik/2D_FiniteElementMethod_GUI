function [X_1, Y_1, stress_xx, stress_xy, stress_yx,stress_yy...
    ,strain_xx,strain_xy,strain_yx,strain_yy, disp_x,disp_y ] = post_process(NL,EL,ENL,scale)
%Plotting and illustration of our results



GPE = getGlobalgpe;

PD = 2;
NPE = size(EL,2);
NoN = size(ENL,1);
NoE = size(EL,1);

    

[disp,stress,strain]= element_pp(NL , EL , ENL );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   PLOTTING    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if NPE == 4
    
    if GPE == 4

         for i = 1:NoE
            %nl = EL(i, 1, NPE);
            nl = EL(i, 1:NPE);
            
            for j = 1:NPE
                X_1(j, i) = ENL(nl(j), 1) +scale*ENL(nl(j), 4*PD+1) ;
                Y_1(j, i) = ENL(nl(j), 2) + scale*ENL(nl(j), 4*PD+2) ;
            end
  %%%%%%%%%%%%%%%%%%%%%%% Stresses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for j = 1:NPE
                val = stress(i, :, 1, 1);
                stress_xx(j,i) = val(1, j);
            end
            
            for j = 1:NPE
                val = stress(i, :, 1, 2);
                stress_xy(j,i) = val(1, j);
            end
            
            for j = 1:NPE
                val = stress(i, :, 2, 1);
                stress_yx(j,i) = val(1, j);
            end
            
            for j = 1:NPE
                val = stress(i, :, 2, 2);
                stress_yy(j,i) = val(1, j);
            end
      %%%%%%%%%%%%%%%%%%%%%%% Displacement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
            for j = 1:NPE
                disp_y(j, i) = abs(ENL(nl(j), 4*PD + 2));
                disp_x(j, i) = abs(ENL(nl(j), 4*PD + 1));
            end
     %%%%%%%%%%%%%%%%%%%%%%% Strains %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for j = 1:NPE
                val2 = strain(i, :, 1, 1);
                strain_xx(j,i) = val2(1, j);
            end
            
            for j = 1:NPE
                val2 = strain(i, :, 1, 2);
                strain_xy(j,i) = val2(1, j);
            end
            
            for j = 1:NPE
                val2 = strain(i, :, 2, 1);
                strain_yx(j,i) = val2(1, j);
            end
            
            for j = 1:NPE
                val2 = strain(i, :, 2, 2);
                strain_yy(j,i) = val2(1, j);
            end
            
        end
    
    elseif GPE == 1
    end
              
    
 elseif (NPE == 3) %Or NPE ==6 or 9 or 10  

        
       for i = 1:NoE

            nl = EL(i,1:NPE); %identify all the nodes connected to that element

            for j=1:NPE  %getting the coordinates of each node
                X_1(j,i) = ENL(nl(j),1) ;
                Y_1(j,i) = ENL(nl(j),2) ;

            end

            for j=1:NPE
                val = stress(i,:,1,1);
                stress_xx(j,i)= val;
            end

            for j = 1:NPE
                val = stress(i, :, 1, 2);
                stress_xy(j,i) = val;
            end

            for j = 1:NPE
                val = stress(i, :, 2, 1);
                stress_yx(j,i) = val;
            end
            
            for j = 1:NPE
                val = stress(i, :, 2, 2);
                stress_yy(j,i) = val;
            end          
            for j = 1:NPE
                disp_y(j, i) = abs(ENL(nl(j), 4*PD + 2));
            end  
            for j = 1:NPE
                disp_x(j, i) = abs(ENL(nl(j), 4*PD + 1));
            end
            
            
            for j = 1:NPE
                val = strain(i, :, 1, 1);
                strain_xx(j,i) = val;
            end
            
            for j = 1:NPE
                val = strain(i, :, 1, 2);
                strain_xy(j,i) = val;
            end
            
            for j = 1:NPE
                val = strain(i, :, 2, 1);
                strain_yx(j,i) = val;
            end
            
            for j = 1:NPE
                val = strain(i, :, 2, 2);
                strain_yy(j,i) = val;
            end

       end

        
  
end
end
    
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [disp, stress, strain] = element_pp(NL, EL, ENL)

GPE = getGlobalgpe;
V = getGlobalV;
E = getGlobalE;
PD = size(NL,2);
NPE = size(EL,2);
NoN = size(ENL,1);
NoE = size(EL,1);


disp   = zeros(NoE, GPE, PD, PD);
stress = zeros(NoE, GPE, PD, PD);



for e=1:NoE %looping over elements

    nl = EL(e,1:NPE); %obtaining all 'NPE' nodes from EL connected to element 'e'

    coor = zeros(PD,NPE); %coordinates of each element
    for i = 1:NPE
        for j = 1:PD
            coor(j,i) = NL(nl(i),j);  
        end
    end
    
    u = zeros(PD, NoE);
    for i = 1:NPE
        for j = 1:PD
            u(j,i) = ENL(nl(i), 4*PD + j);
        end
    end
    for i = 1:NPE
        for j = 1:PD
            disp(e,i,j,1) = ENL(nl(i), 4*PD + j);
        end
    end

    for gp = 1:GPE
        epsilon = zeros(PD, PD);
        for i = 1:NPE
            J = zeros(PD, PD);
            grad = zeros(PD, NPE);
        end
    end

    for gp = 1:GPE %loop over Gauss Points

        epsilon = zeros(PD,PD);  

        for i=1:NPE

                J = zeros(PD,PD);  %JACOBIAN

                grad = zeros(PD,PD);

                [xi,eta,alpha] = GaussPoint(NPE,GPE,gp);

                grad_nat = grad_N_nat( NPE, xi , eta);

                J = coor * grad_nat';

                grad = inv(J)' * grad_nat;

                   
                epsilon = epsilon + 1/2 * (dyad(u(:, i), grad(:, i)) + dyad(grad(:, i), u(:, i)));
        end
                sigma = zeros(PD, PD);
        for a = 1:PD
            for b = 1:PD
                for c = 1:PD
                    for d = 1:PD
                        sigma(a, b) = sigma(a, b) + constitutive(a, b, c, d,E,V) * epsilon(c, d);
                    end
                end
            end
        end
        
        for a = 1:PD
            for b = 1:PD
                strain(e,gp,a,b) = epsilon(a,b);
                stress(e,gp,a,b) = sigma(a,b);
            end
        end

        end  

 end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        result(2,7) = (1/2)*(1+xi)*(1-xi);
        result(2,8) = -eta*(1-xi);
        
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ xi , eta, alpha ] = GaussPoint( NPE , GPE , gp)
    GPE = getGlobalgpe;
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
        if (GPE == 9)
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
        elseif (GPE == 1) %if 1 gauss points
            xi = 0 ; eta = 0; alpha = 4;
        end
    end
end

function A = dyad(u, v)



PD = 2;
A = zeros(PD, PD);
for i = 1:PD
    for j = 1:PD
        A(i, j) = u(i) * v(j);
    end
end


end

function C = constitutive(i, j, k, l,E,V)

C = (E/(2*(1+V))) * ((i==l)*(j==k) + (i==k)*(j==l)) + E*V/(1-V^2) * (i==j)*(k==l);
end


function GPE = getGlobalgpe
global gpe
GPE = gpe;
end


function v = getGlobalV
global V
v = V;
end

function e = getGlobalE
global E
e = E;
end

