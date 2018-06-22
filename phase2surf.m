function [ Zsurf ] = phase2surf( p,z_light,z_sensor,x1,x2,y1,y2,h,A1,A2,B1,B2 )
%==========================================================================
%function used data collected when parallel light rays shoot down onto a
%to reconstruct an unknown flat surface
%-----------------------------Inputs---------------------------------------
%p= nxm matrix that contains the p(x,y) values of the phase function
%z_light= height of the light source
%z_sensor= height of the light capturing sensor
%[x1,x2]x[y1,y2] = domain of the light source on a coordinate plane
%h= step size of the mesh
%[A1,A2] x [B1,B2] = domain of the light when it reflects onto the sensor
%-----------------------------Output---------------------------------------
%Zsurf= nxm matrix that approximates the unknown surface
%---------------------------Assumptions------------------------------------
%1. The reflected light rays that hits the sensor will hit it in such a
%   way that the domain that is captured is will be a rectangle.
%2. The light ray shining down have the direction w=[0,0,-1].
%3. The elements for the gradient of p doesn't no equal zero, if so there
%   will be an instance of # divided by zero.
%==========================================================================
[n,m] = size(p);
Zsurf=zeros(n,m);
% px = zeros(n,m);
% py = zeros(n,m);
% 
% for i = 2:n-1
%     for j = 2:m-1
%         px(i,j)= (p(i+1,j) - p(i-1,j))/2*h;
%         py(i,j)= (p(i,j+1) - p(i,j-1))/2*h;
%     end
% end
[px,py]=gradient(p);

%make the mesh for the domain on the surface
X=x1:h:x2;
Y=y1:h:y2;
sx=A1:(A2-A1)/(n-1):A2;
sy=B1:(B2-B1)/(m-1):B2;

%==========get the height of the surface for the boundaries================
% -------approximate the z=f(x,y) for the upper and lower boundaries-------
for j=1:m
    tBy=(sy(j)-Y(j))/(py(1,j));
    tBx=(sx(1)-X(1))/px(1,j);
%     tB=(tBx+tBy)/2
    Zsurf(1,j)=z_sensor-max(tBx,tBy);
    
    tTy=(sy(j)-Y(j))/(py(n,j));
    tTx=(sx(n)-X(n))/px(n,j);
%     tT=(tTy+tTx)/2
    Zsurf(n,j)=z_sensor-max(tTx,tTy);
end

%------approximate the z=f(x,y) for the left and right boundaries----------
for i=2:n-1
    tLx=(sx(i)-X(i))/(px(i,1));
    tLy=(sy(1)-Y(1))/py(i,1);
    Zsurf(i,1)=z_sensor-max(tLx,tLy);
    
    tRx=(sx(i)-X(i))/(px(i,m));
    tRy=(sy(m)-Y(m))/py(i,m);
    Zsurf(i,m)=z_sensor-max(tRx,tRy);
end


%=====================choose favorite k value==============================
%get the k values for the corners of the domain and take the min
k1=z_light+(z_sensor-Zsurf(1,1))+ norm(([X(1),Y(1),Zsurf(1,1)]-[A1,B1,z_sensor]));
k2=z_light+(z_sensor-Zsurf(1,m))+ norm(([X(1),Y(m),Zsurf(1,m)]-[A1,B2,z_sensor]));
k3=z_light+(z_sensor-Zsurf(n,1))+ norm(([X(n),Y(1),Zsurf(n,1)]-[A2,B1,z_sensor]));
k4=z_light+(z_sensor-Zsurf(n,m))+ norm(([X(n),Y(m),Zsurf(n,m)]-[A2,B2,z_sensor]));
kk=[k1,k2,k3,k4];
k=min(kk);

%-------------------get the inner height values----------------------------
w=[0,0,-1];
nw=norm(w);
for i=1:n
    for j=1:m
        gphi=[px(i,j),py(i,j),1];
        uw=abs(dot(gphi,w)/(norm(gphi)*nw));
        h1=z_light-z_sensor;
        h3=(k-h1)/(uw+1);
        h2=abs(h3*uw);
        Zsurf(i,j)= z_sensor-h2;
%         h3=z_light-z_sensor-p(i,j);
%         h1=h3-k-uw;
%         h2=h1*uw;
%         Zsurf(i,j)=z_sensor-h2;
    end
end
Zsurf=Zsurf+p;
end

