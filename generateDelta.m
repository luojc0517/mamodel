function [ delta ] = generateDelta( X,beta )
%GENERATEDELTA X,beta
%  得到delta三角形
 mu=X*beta;
 delta=diag(mu);
end

