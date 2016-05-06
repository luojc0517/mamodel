function [ gamma_new ] = updateGamma1( Y,mu,Z,L,m,t,var,D,gamma_old)
%[ gamma_new ] = updateGamma1( Y,mu,L,m,t,var,D,gamma_old)
%   迭代更新gamma,主意mu不能是初始值，而是应该使用beta迭代后的mu
gamma0=[-0.3,0.3,0,0,0,0,0]';
n=length(Y);
R=cell(1,n);
for i=1:n
    Ri=Y{i}-mu{i};
    R{i}=Ri;
end
E=cell(1,n);
E{1}=R{1}';
for i=2:n
    Li=L{i};
    for j=1:m{i}
        sum=0;
        for k=1:j-1
            sum=sum+Li(j,k)*E{i}(k);
        end
        E{i}(j)=R{i}(j)-sum;
    end
end
dev_e_gamma=cell(1,n);
for i=1:n
    dev_e_gamma{i}=[];
    dev_e_gamma{i}(:,1)=zeros(length(gamma0),1);
    for j=2:m{i}
        sum=0;
        for k=1:j-1
            [ zijk,lijk ] = generateLijk( t{i},gamma0,j,k );
            sum=sum+E{i}(k)*zijk'+L{i}(j,k)*dev_e_gamma{i}(:,k);
        end
        dev_e_gamma{i}(:,j)=-sum;
    end
end
U_gamma=0;
for i=1:n
    U_gamma=U_gamma+dev_e_gamma{i}*inv(D{i})*E{i}';
end
sum_Ggamma=0;
for i=1:n
    [ G_gamma ] = generateGgamma( L{i},Z{i},var{i},m{i} );
    sum_Ggamma=sum_Ggamma+G_gamma;
end
Ggamma_sum=0;
for i=1:n
    [ G_gamma ] = generateGgamma( L{i},Z{i},var{i},m{i} );
    Ggamma_sum=Ggamma_sum+G_gamma;
end

gamma_new=gamma_old-inv(Ggamma_sum)*U_gamma;



end

