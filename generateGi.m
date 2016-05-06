A=inv(L);
G_gamma=0;
for j=2:m
    sum2=0;
   for k=1:j-1
       sum2=sum2+A(j,k)*Z{k};
   end
   part2=Z{j}+sum2;
   G_gamma=G_gamma+(1./var')*part2*D*part2';
end