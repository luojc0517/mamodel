function [ newBeta,mu ] = updateBeta(beta,lamda,gamma,X,Y,H,t)
%UPDATEBETA beta,lamda,gamma,X,Y,H,t
%   更新beta
n=length(X);%样本个数
delta=cell(1,n);
sigma=cell(1,n);
mu=cell(1,n);
sum=zeros(10,1);
for i=1:n
    delta{i}=generateDelta(X{i},beta);
    sigma{i}=generateSigma( lamda,gamma,t{i},H{i});
    mu{i}=X{i}*beta;
    sum=sum+X{i}'*delta{i}*(sigma{i}\(Y{i}-mu{i}));
end
IBeta=generateIBeta( X,delta,sigma);
newBeta=beta+IBeta\sum;
end

