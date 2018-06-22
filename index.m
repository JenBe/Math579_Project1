function [ num ] = index( n, i,j )
%Takes the i,j index of the an element and converts it to a number

%Inputs:
%   n - length of the row
%   i,j - index value

%Output:
%   num - represents that i,j is the kth element of the matrix

num=(i-1)*n+j;
end

