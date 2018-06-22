function [ reflection_x,reflection_y,reflection_z ] = getReflections( surface,surface_normal,incident_direction,x1,x2,y1,y2,h, t )
%--------------------------PURPOSE-----------------------------------------
%gets the reflected values
%--------------------------INPUTS------------------------------------------
%surface = function handle of the surface
%surface_normal = function handle of the gradient
%incident_direction = where the light rays are shooting from
%[x1,x2]x[y1,y2] = domain
%h= step size
%t = "time"
%-------------------------OUTPUT-------------------------------------------
%reflected x,y,z coordinates
%--------------------------------------------------------------------------
x=x1:h:x2; y=y1:h:y2;
n=length(x); m=length(y);

reflection_x=zeros(n,m);
reflection_y=reflection_x;
reflection_z=reflection_x;

for i= 1:n
    for j=1:m
        z=surface(x(i),y(j));
        temp_normal=surface_normal(x(i),y(j),z);
        gradientphi=temp_normal-incident_direction;
        
        reflection=[x(i),y(j),z] + gradientphi*t;
        reflection_x(i,j)=reflection(1);
        reflection_y(i,j)=reflection(2);
        reflection_z(i,j)=reflection(3);
    end
end

end

