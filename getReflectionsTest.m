%%spherical surface
% r=2;
% surface=@(x,y) sqrt((r^2)-(x^2)-(y^2));
% surface_normal=@(x,y,z) [2*x,2*y,2*z];
% incident_direction=[0,0,-1];
% 
% x1=-1; x2=1; y1=-1; y2=1;
% h=0.2;
% x=x1:h:x2; y=y1:h:y2;
% [X,Y]=meshgrid(x,y);
% n=length(x);m=length(y);
% Zexact=zeros(n,m);
% for i=1:n
%     for j=1:m
%         Zexact(i,j)=surface(x(i),y(j));
%     end
% end
% 
% [xt1,yt1,zt1]= getReflections( surface,surface_normal,incident_direction,x1,x2,y1,y2,h, 1 );
% [xt2,yt2,zt2]= getReflections( surface,surface_normal,incident_direction,x1,x2,y1,y2,h, 2 );
% [xt3,yt3,zt3]= getReflections( surface,surface_normal,incident_direction,x1,x2,y1,y2,h, 3 );
% 
% p=zt2-2;
% z_light=20;
% z_sensor=3;
% A1=-7;A2=7;B1=-7;B2=7;
% Zsurf = phase2surf( p,z_light,z_sensor,x1,x2,y1,y2,h,A1,A2,B1,B2 );
% subplot(1,3,1)
% surf(X,Y,Zsurf)
% title('Approximated Sphere')
% subplot(1,3,2)
% surf(X,Y,Zexact)
% title('Exact Sphere')
% 
% subplot(1,3,3)
% mesh(xt1,yt1,zt1)
% title('Reflected Spherical Wavefronts')
% hold on
% mesh(xt2,yt2,zt2)
% mesh(xt3,yt3,zt3)
% hold off
% 
% 
% % ellipsoid surface
% A=3; B=4; C=2;
% surface=@(x,y) sqrt((C^2)*(1-((x^2)/(A^2))-((y^2)/(B^2))));
% surface_normal=@(x,y,z) [2*x/(A^2),2*y/(B^2),2*z/(C^2)];
% incident_direction=[0,0,-1];
% 
% x1=-1; x2=1; y1=-1; y2=1;
% h=0.2;
% x=x1:h:x2; y=y1:h:y2;
% [X,Y]=meshgrid(x,y);
% 
% n=length(x);m=length(y);
% Zexact=zeros(n,m);
% for i=1:n
%     for j=1:m
%         Zexact(i,j)=surface(x(i),y(j));
%     end
% end
% 
% [xt1,yt1,zt1]= getReflections( surface,surface_normal,incident_direction,x1,x2,y1,y2,h, 1 );
% [xt2,yt2,zt2]= getReflections( surface,surface_normal,incident_direction,x1,x2,y1,y2,h, 2 );
% [xt3,yt3,zt3]= getReflections( surface,surface_normal,incident_direction,x1,x2,y1,y2,h, 3 );
% 
% p=zt3-zt2;
% z_light=20; z_sensor=3;
% A1=-5/3; A2=5/3; B1=-1.3750; B2=1.3750;
% Zsurf = phase2surf( p,z_light,z_sensor,x1,x2,y1,y2,h,A1,A2,B1,B2 );
% 
% subplot(1,3,1)
% surf(X,Y,Zsurf)
% title('Approximated Ellipsoid')
% subplot(1,3,2)
% surf(X,Y,Zexact)
% title('Exact Ellipsoid')
% hold on
% subplot(1,3,3)
% mesh(xt1,yt1,zt1)
% title('Reflected Ellipsoid Wavefronts')
% hold on
% mesh(xt2,yt2,zt2)
% hold on
% mesh(xt3,yt3,zt3)
% hold off

%% plane surface
surface=@(x,y) .2*x+.2*y+1;
surface_normal=@(x,y,z) [.2,+.2,1];

x1=-1; x2=1; y1=-1; y2=1;
h=0.1;
x=x1:h:x2; y=y1:h:y2;
[X,Y]=meshgrid(x,y);

n=length(x);m=length(y);
Zexact=zeros(n,m);
for i=1:n
    for j=1:m
        Zexact(i,j)=surface(x(i),y(j));
    end
end

% Getting the reflected srufaces from known surface
[xt1,yt1,zt1]= getReflections( surface,surface_normal,incident_direction,x1,x2,y1,y2,h, 1 );
[xt2,yt2,zt2]= getReflections( surface,surface_normal,incident_direction,x1,x2,y1,y2,h, 2 );
[xt3,yt3,zt3]= getReflections( surface,surface_normal,incident_direction,x1,x2,y1,y2,h, 3 );


% Getting back our unknown surface
p=zt3+1;
z_light=10; z_sensor=3;
A1=-0.25; A2=1.75; B1=-2.5; B2=-0.5;
Zsurf = phase2surf( p,z_light,z_sensor,x1,x2,y1,y2,h,A1,A2,B1,B2 );

subplot(1,3,1)
surf(X,Y,Zsurf)
title('Approximated Uneven Plane')
subplot(1,3,2)
surf(X,Y,Zexact)
title('Exact Uneven Plane')
hold on
subplot(1,3,3)
mesh(xt1,yt1,zt1)
title('Reflected Uneven Plane Wavefronts')
hold on
mesh(xt2,yt2,zt2)
hold on
mesh(xt3,yt3,zt3)
hold off