function [X10,X10_pre,HCS1,MeanAE,MeanRE,AIC,NASH]=GM13Prediction(runoff,temperature,precipitation,yearb,yeare,t_pre,k)
data1=runoff';
data2=temperature';
data3=precipitation';
%-----��ֵ��׼��------
X101=data1(1);data1=data1./X101;
X201=data2(1);data2=data2./X201;
X301=data3(1);data3=data3./X301;
%-----------------
n=length(data1);  
X10=data1;
X20=data2;
X30=data3;
%������ʼ��
X11=zeros(1,n);X21=zeros(1,n);X31=zeros(1,n);
M1=zeros(1,n-1);M2=zeros(1,n-1);M3=zeros(1,n-1);
Y1=zeros(1,n-1);Y2=zeros(1,n-1);Y3=zeros(1,n-1);
X31_pre=zeros(1,n+t_pre);X30_pre=zeros(1,n+t_pre);
X21_pre=zeros(1,n+t_pre);X20_pre=zeros(1,n+t_pre);
X11_pre=zeros(1,n+t_pre);X10_pre=zeros(1,n+t_pre);
rawX10_pre=zeros(1,n+t_pre);
%�ۼ����ɴ���
X11(1)=X10(1);
X21(1)=X20(1);
X31(1)=X30(1);
for i=2:n  
   X11(i)=X11(i-1)+X10(i); 
   X21(i)=X21(i-1)+X20(i);
   X31(i)=X31(i-1)+X30(i);
end 
%������ھ�ֵ��������
for i=1:n-1 
   M1(i)=(0.5*(X11(i)+X11(i+1)));
   M2(i)=(0.5*(X21(i)+X21(i+1)));
   M3(i)=(0.5*(X31(i)+X31(i+1)));
end 
 B1=zeros(n-1,3); 
for i=1:(n-1) 
    B1(i,1)=-M1(i);    
    B1(i,2)=X21(i+1); 
    B1(i,3)=X31(i+1);
end
B2=zeros(n-1,2); 
for i=1:(n-1) 
    B2(i,1)=-M2(i);  
    B2(i,2)=X31(i+1); 
end
B3=zeros(n-1,2); 
for i=1:(n-1) 
    B3(i,1)=-M3(i);   
    B3(i,2)=1; 
end
%���쳣��������Y
for i=2:n                          
    Y1(i-1)=X10(i); 
    Y2(i-1)=X20(i);
    Y3(i-1)=X30(i);
end 
HCS1=inv(B1'*B1)*B1'*Y1';               %����С���˷���Ҳ���HCS1 
H1=HCS1';                            %H1=[a,b2,b3]
HCS2=inv(B2'*B2)*B2'*Y2';               %����С���˷���Ҳ���HCS2 
H2=HCS2';                            %H2=[a,b3]
HCS3=inv(B3'*B3)*B3'*Y3';               %����С���˷���Ҳ���HCS3 
H3=HCS3';                            %H3=[a,b]
%�����X3���ۼ�����
for i=1:n+t_pre                         
X31_pre(i)=(X30(1)-H3(2)/H3(1))*exp(-1*H3(1)*(i-1))+H3(2)/H3(1); 
end 
for i=2:n+t_pre                      
       X30_pre(i)=X31_pre(i)-X31_pre(i-1);
end
X30_pre(1)=X30(1);

%�Բ�����alpha��beta�任
H2=H2./(1+0.5*H2(1));
%��ԭ�����X2��Ԥ��ֵ
X20_pre(1)=X20(1);
for i=2:n                     
       X20_pre(i)=H2(2).*X30(i)+(1-H2(1)).*X20(i-1);
end
X21_pre(n)=X21(n);
for i=n+1:n+t_pre
    X20_pre(i)=H2(2).*X30_pre(i)+(1-H2(1)).*X20_pre(i-1);
    X21_pre(i)=X20_pre(i)+X21_pre(i-1);
end
%�Բ�����alpha��beta�任
H1=H1./(1+0.5*H1(1));
%��ԭ�����X1��Ԥ��ֵ
X10_pre(1)=X10(1);
for i=2:n                     
       X10_pre(i)=H1(2).*X20(i)+H1(3).*X30(i)+(1-H1(1)).*X10(i-1);
end
X11_pre(n)=X11(n);
for i=n+1:n+t_pre
    X10_pre(i)=H1(2).*X20_pre(i)+H1(3).*X30_pre(i)+(1-H1(1)).*X10_pre(i-1);
    X11_pre(i)=X10_pre(i)+X11_pre(i-1);
end
%------��ֵ���׼��--------
X10_pre=X10_pre.*X101;
X20_pre=X20_pre.*X201;
X30_pre=X30_pre.*X301;
X10=X10.*X101;
X20=X20.*X201;
X30=X30.*X301;
%%
%��ͼ
xraw=yearb:yeare;
yraw=X10(1:n);
xpre=yearb+1:yeare+t_pre;
ypre=X10_pre(2:n+t_pre);
plot(xraw,yraw,'-*r'),xlabel('���'),ylabel('����'),title('ԭʼ������Ԥ��ֵ�Ա�');
hold on;
plot(xpre,ypre,'-ob');
% legend('����','����','��ˮ','����Ԥ��ֵ')
legend('����','����Ԥ��ֵ')
%%
%����
%�������
AE=X10(1:n)-X10_pre(1:n);%AE=AbsoluteError
%ƽ���������
MeanAE=mean(abs(AE));
%������
RE=AE./X10(1:n);%RE=RelativeError
%ƽ��������
MeanRE=mean(abs(RE));
%�����c
AveX10=mean(X10);%Ave=Average
VarX10= sum((X10-AveX10).^2)/n;%Variance
AveAE=mean(AE(2:n));
VarAE=sum((AE(2:n)-AveAE).^2)/(n-1);
c=VarAE/VarX10;
%С������p
e=abs(AE(2:n)-AveAE);
S=0.6745*sqrt(VarX10);
p=sum(e<S)/(n-1);
%NASHϵ��
NASH=1-sum((X10(2:n)-X10_pre(2:n)).*(X10(2:n)-X10_pre(2:n)))/sum((X10(2:n)-mean(X10(2:n))).*(X10(2:n)-mean(X10(2:n))));
%AIC ��С��Ϣ׼��
% k=4;%ע���޸�
SSE=sum((X10(2:n)-X10_pre(2:n)).*(X10(2:n)-X10_pre(2:n)));%�в�ƽ����
l=-(n/2)*log(2*pi)-(n/2)*log(SSE/n)-n/2;
AIC=(-2*l+2*k)/n; %�׵��Ե�<<���ݷ�����EVIEWSӦ��>>














