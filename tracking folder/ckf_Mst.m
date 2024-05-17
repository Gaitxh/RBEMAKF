function z=ckf_Mst(x,nz,H)

n=size(x,2);

for i=1:n
    z(:,i)=MstEq(x(:,i),nz,H);
end
