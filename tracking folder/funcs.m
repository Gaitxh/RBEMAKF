function [fx]=funcs(gama,n,sigma,v,flag)

%%%%%%%��һ��ѡ��
if flag==1
    fx = sigma*sigma*exp(0.5*(n-gama)/sigma^2);
end

%%%%%%%�ڶ���ѡ��
if flag==2
    fx = -0.5*(v+n)*log(n + gama/v);
end

%%%%%%%������ѡ��
if flag==3
    fx = -sqrt((v+n)*(v+gama));
end


end

