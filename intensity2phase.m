%takes the intensity and give a phase function with a the desired z
%inputs same as TIENeumann, z
%output: phi = phase function

function [ phi ] = intensity2phase( I,f,h,sigma,z )

%use TIENeumann to get p(x,y)
[p,lambda]=TIENeumann(I,f,h,sigma);

%phase function = z + p(x,y)
phi=z+p;

end

