function UM=Build_Mat(u,i,jj,N)
mp=size(u,2)-1;
UM=zeros(mp*N,jj);
cont=1;
for k=i:i+jj-1;
    dummyV1=zeros(mp*N,1);
    for j=k:N+k-1
        dummyV2=u(j,2:end)';
        dummyV1((j-k)*mp+1:(j-k+1)*mp)=dummyV2;
    end
    UM(:,cont)=dummyV1;
    cont=cont+1;
end