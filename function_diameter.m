function [d32, d50] = function_diameter(A)
%GRANULOMETRY gives the granulometric analysis of the particulate system.
%
%   Z = GRANULOMETRY(H1,CLASES) returns the granulometric report of h1 in a
%   given number of clases "CLASES". H1 is a vector.
%   Abhinandan Singh 16.01.2018
smd=zeros(1,2);
[uvals, ~, uidx] = unique(A);
output = [uvals, accumarray(uidx, 1)];
for i=1:1:size(output)
    smd(i,1)=output(i,2)*(output(i,1).^3);
    smd(i,2)=output(i,2)*(output(i,1).^2);
end
d32=sum(smd(:,1))/sum(smd(:,2));  
dmax=max(A);
dmin=min(A);
d50=dmin+(0.32*(dmax-dmin));
end
