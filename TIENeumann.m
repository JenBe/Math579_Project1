function [p, lambda] = TIENeumann( I, f, h, sigma )
%TIENeumann Function that solves the Transport of Intensity Equation with
%homogeneous Nemann boundary condition.

%Inputs:
%   I - 2D array that stores the intensity data on D = D U dD
%   f - 2D array that has the values f(x,y) on D
%   h - number that represents mesh size
%   sigma - number that represents the sum of p(x,y) over D

%Outputs:
%   p - 2D array that stores the p(x,y) values on D
%   lambda - number that measure the validity of the compatibility
%           condition (double integral of f dxdy over D = 0)

[n,m] = size(f);
F = zeros(m*n,1);
A= zeros(m*n,m*n);

%--------------------------Build A Matrix----------------------------------
for i=2:n-1
    for j=2:m-1
        ind=index(m,i,j);
        ind1=index(m,i-1,j);
        ind2=index(m,i,j-1);
        ind4=index(m,i,j+1);
        ind5=index(m,i+1,j);
        
       A(ind,ind1)=I(i,j)+I(i-1,j);
       A(ind,ind2)=I(i,j)+I(i,j-1);
       A(ind,ind)=-(I(i,j+1)+I(i,j-1)+I(i+1,j)+I(i-1,j)+4*I(i,j));
       A(ind,ind4)=I(i,j+1)+I(i,j);
       A(ind,ind5)=I(i+1,j)+I(i,j);
    end
end

%---------------insert values for the edges ----------------------------
% Neumann conditions for edges where top and bottom edges are evaluated
% with -I*partial(P)/partial(y) from corner to corner (1,1)....(1,m) for 
% bottom and I*partial(P)/partial(y) (n,1)....(n,m) for top. The side edges
% are evaluated at -I*partial(P)/partial(x) for left edge from
% (2,1)...(n-1,1) and I*partial(P)/partial(x) from (2,m)...(n-1,m) for
% right edge.
%********** Upper/Lower Boundary ****** (+/-)I*partial(P)/partial(y)*****
for j = 1:m
    % Upper Edge ----------
    indUpper=index(m,n,j);
    ind1 = index(m,n,j);
    ind2 = index(m,n-1,j);
    ind3 = index(m,n-2,j);
    coef = I(n,j);
    A(indUpper,ind1) = -3*coef;
    A(indUpper,ind2) = 4*coef;
    A(indUpper,ind3) = -1*coef;
    
    % Lower Edge ------------
    indLower = index(m,1,j);
    index1 = index(m,1,j); 
    index2 = index(m,2,j);
    index3 = index(m,3,j);
    num = I(1,j);
    
    A(indLower,index1) = -3*num;
    A(indLower,index2) = 4*num;
    A(indLower,index3) = -num;
end
%****Left/Right Boundaries *******(+/-)I*partial(P)/partial(x)*********
% These points do not include the corners since they are filled in as part
% of the upper and lower boundaries.
for i = 2:(n-1)
    % Left boundary ----------------------
    indLeft=index(m,i,1);
    indL1=index(m,i,1);
    indL2=index(m,i,2);
    indL3=index(m,i,3);
    xCoeff = I(i,1);
    A(indLeft,indL1) = (-3*xCoeff);
    A(indLeft,indL2) = (4*xCoeff);
    A(indLeft,indL3)= -xCoeff;
    
    % Right Boundary ----------------------
    indRight=index(m,i,m);
    indR1 = index(m,i,m);
    indR2 = index(m,i,m-1);
    indR3 = index(m,i,m-2);
    rxCoeff = I(i,m);
    A(indRight,indR1)= (-3*rxCoeff);
    A(indRight,indR2) = 4*rxCoeff;
    A(indRight,indR3) = -rxCoeff;
    
end

%Build A_bar
OnesV=ones(m*n,1);
A_bar=[A,OnesV; OnesV',0];


%---------------------------Build F Vector--------------------------------
for i = 1:n
    for j = 1:m
        ind = index(m,i,j);
        F(ind) = (2*h^2).*f(i,j);
    end
end

%insert values for the edges
% Upper/Lower edges ---------
for k = 1:m
    F(index(m,1,k))=0;
    F(index(m,n,k))=0;
end

% Left/Right edges -----------
for q = 1:n
    F(index(m,q,1))=0;
    F(index(m,q,m))=0;
end

%create F_bar vector
F_bar=[F;sigma];

%-----------------------------Get Outputs---------------------------------
%Solve for P_bar
P_bar=A_bar\F_bar;

%Make p matrix
p=zeros(n,m);
for i=1:n
    for j=1:m
        ind=index(m,i,j);
        p(i,j)=P_bar(ind);
    end
end

%Get lambda value
lambda=P_bar(m*n+1);
end
