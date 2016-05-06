function [ IBeta ] = generateIBeta( X,delta,sigma)
%GENERATEIBETA X,delta,sigma,n
%   µÃµ½IBeta
IBeta=zeros(10,10);
n=length(X);
for i=1:n
IBeta=IBeta+X{i}'*delta{i}*(sigma{i}\delta{i}*X{i});
end
end

