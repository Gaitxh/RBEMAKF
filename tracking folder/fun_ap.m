function [lamda1,lamda2]=fun_ap(gama1,gama2,sigma,v,flag,nx,nz)

%%%%%%%��һ��ѡ��
if flag==1
    lamda1=exp(0.5*(nx-gama1)/sigma^2);
    lamda2=exp(0.5*(nz-gama2)/sigma^2);
end

%%%%%%%�ڶ���ѡ��
if flag==2
    lamda1=(v+nx)/(v+gama1);
    lamda2=(v+nz)/(v+gama2);
end

%%%%%%%������ѡ��
if flag==3
    lamda1=sqrt((v+nx)/(v+gama1));
    lamda2=sqrt((v+nz)/(v+gama2));
end

%%%%%%%��ֹ����
if lamda1<1e-8
    lamda1=lamda1+1e-8;
end
%%%%%%%��ֹ����
if lamda2<1e-8
    lamda2=lamda2+1e-8;
end

end

