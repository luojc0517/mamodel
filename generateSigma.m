function [ Sigma,L,var ] = generateSigma( lamda,gamma,t,H )
%GENERATESIGMA lamda,gamma,t,H
%   ͨ��lamda��gamma�õ�Sigma
 m=length(t);
 L=eye(m)*1;
 Z=[];
 for j=2:m
     for k=1:j-1
         [zijk,lijk]=generateLijk(t,gamma,j,k);
         L(j,k)=lijk;
     end
 end
 
 
 logVar=H*lamda;
 var=exp(logVar);
 D=diag(var);
 Sigma=L*D*L';
end

