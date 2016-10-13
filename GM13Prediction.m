function [X10,X10_pre,HCS1,MeanAE,MeanRE,AIC,NASH]=GM13Prediction(runoff,temperature,precipitation,yearb,yeare,t_pre,k)
data1=runoff';
data2=temperature';
data3=precipitation';
%-----初值标准化------
X101=data1(1);data1=data1./X101;
X201=data2(1);data2=data2./X201;
X301=data3(1);data3=data3./X301;
%-----------------
n=length(data1);  
X10=data1;
X20=data2;
X30=data3;
%变量初始化
X11=zeros(1,n);X21=zeros(1,n);X31=zeros(1,n);
M1=zeros(1,n-1);M2=zeros(1,n-1);M3=zeros(1,n-1);
Y1=zeros(1,n-1);Y2=zeros(1,n-1);Y3=zeros(1,n-1);
X31_pre=zeros(1,n+t_pre);X30_pre=zeros(1,n+t_pre);
X21_pre=zeros(1,n+t_pre);X20_pre=zeros(1,n+t_pre);
X11_pre=zeros(1,n+t_pre);X10_pre=zeros(1,n+t_pre);
rawX10_pre=zeros(1,n+t_pre);
%累加生成处理
X11(1)=X10(1);
X21(1)=X20(1);
X31(1)=X30(1);
for i=2:n  
   X11(i)=X11(i-1)+X10(i); 
   X21(i)=X21(i-1)+X20(i);
   X31(i)=X31(i-1)+X30(i);
end 
%构造紧邻均值生成序列
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
%构造常数项向量Y
for i=2:n                          
    Y1(i-1)=X10(i); 
    Y2(i-1)=X20(i);
    Y3(i-1)=X30(i);
end 
HCS1=inv(B1'*B1)*B1'*Y1';               %用最小二乘法求灰参数HCS1 
H1=HCS1';                            %H1=[a,b2,b3]
HCS2=inv(B2'*B2)*B2'*Y2';               %用最小二乘法求灰参数HCS2 
H2=HCS2';                            %H2=[a,b3]
HCS3=inv(B3'*B3)*B3'*Y3';               %用最小二乘法求灰参数HCS3 
H3=HCS3';                            %H3=[a,b]
%计算出X3的累加序列
for i=1:n+t_pre                         
X31_pre(i)=(X30(1)-H3(2)/H3(1))*exp(-1*H3(1)*(i-1))+H3(2)/H3(1); 
end 
for i=2:n+t_pre                      
       X30_pre(i)=X31_pre(i)-X31_pre(i-1);
end
X30_pre(1)=X30(1);

%对参数作alpha，beta变换
H2=H2./(1+0.5*H2(1));
%还原计算出X2的预测值
X20_pre(1)=X20(1);
for i=2:n                     
       X20_pre(i)=H2(2).*X30(i)+(1-H2(1)).*X20(i-1);
end
X21_pre(n)=X21(n);
for i=n+1:n+t_pre
    X20_pre(i)=H2(2).*X30_pre(i)+(1-H2(1)).*X20_pre(i-1);
    X21_pre(i)=X20_pre(i)+X21_pre(i-1);
end
%对参数作alpha，beta变换
H1=H1./(1+0.5*H1(1));
%还原计算出X1的预测值
X10_pre(1)=X10(1);
for i=2:n                     
       X10_pre(i)=H1(2).*X20(i)+H1(3).*X30(i)+(1-H1(1)).*X10(i-1);
end
X11_pre(n)=X11(n);
for i=n+1:n+t_pre
    X10_pre(i)=H1(2).*X20_pre(i)+H1(3).*X30_pre(i)+(1-H1(1)).*X10_pre(i-1);
    X11_pre(i)=X10_pre(i)+X11_pre(i-1);
end
%------初值逆标准化--------
X10_pre=X10_pre.*X101;
X20_pre=X20_pre.*X201;
X30_pre=X30_pre.*X301;
X10=X10.*X101;
X20=X20.*X201;
X30=X30.*X301;
%%
%绘图
xraw=yearb:yeare;
yraw=X10(1:n);
xpre=yearb+1:yeare+t_pre;
ypre=X10_pre(2:n+t_pre);
plot(xraw,yraw,'-*r'),xlabel('年份'),ylabel('径流'),title('原始数据与预测值对比');
hold on;
plot(xpre,ypre,'-ob');
% legend('径流','气温','降水','径流预测值')
legend('径流','径流预测值')
%%
%检验
%绝对误差
AE=X10(1:n)-X10_pre(1:n);%AE=AbsoluteError
%平均绝对误差
MeanAE=mean(abs(AE));
%相对误差
RE=AE./X10(1:n);%RE=RelativeError
%平均相对误差
MeanRE=mean(abs(RE));
%方差比c
AveX10=mean(X10);%Ave=Average
VarX10= sum((X10-AveX10).^2)/n;%Variance
AveAE=mean(AE(2:n));
VarAE=sum((AE(2:n)-AveAE).^2)/(n-1);
c=VarAE/VarX10;
%小误差概率p
e=abs(AE(2:n)-AveAE);
S=0.6745*sqrt(VarX10);
p=sum(e<S)/(n-1);
%NASH系数
NASH=1-sum((X10(2:n)-X10_pre(2:n)).*(X10(2:n)-X10_pre(2:n)))/sum((X10(2:n)-mean(X10(2:n))).*(X10(2:n)-mean(X10(2:n))));
%AIC 最小信息准则
% k=4;%注意修改
SSE=sum((X10(2:n)-X10_pre(2:n)).*(X10(2:n)-X10_pre(2:n)));%残差平方和
l=-(n/2)*log(2*pi)-(n/2)*log(SSE/n)-n/2;
AIC=(-2*l+2*k)/n; %易丹辉的<<数据分析与EVIEWS应用>>














