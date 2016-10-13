%-----------------------------------------------------------------
%GM(1,3)模型
%请直接运行GM13_Main.m即可。
%原理：建立累加序列X1，通过最小二乘法估计灰微分方程的参数，将参数代入白化动态方程还原X0；
%详见华中理工大学出版社1990年出版的邓聚龙《多维灰色规划》第一章GM（1，N)模型；
%数据输入描述：
%runoff.txt为待预测因子runoff纵向排列的时间序列数据；
%temperature.txt为影响因子temperature纵向排列的时间序列数据；
%precipitation.txt为影响因子precipitation纵向排列的时间序列数据；
% YearBegin为开始年份；
% YearEnd为终止年份；
% Years_pre为待预测的年数；
% k为GM13的参数数量，显然等同于GM13中的3；
%数据输出描述：
% X10为待预测因子runoff的实际观测值；
% X10_pre为待预测因子runoff的模拟值；
% HCS1为GM13模型的参数a，b1，b2；分别是runoff，temperature，precipitation的参数。
% MeanAE顾名思义是平均绝对误差
% MeanRE是平均相对误差
% AIC为最小信息准则，是参数数量k和预测残差的函数，越小越好；
% NASH系数为水文学领域描述模型确定性的系数；
% Author: Chong Wang E-mail:chong0314@gmail.com
%------------------------------------------------------------------
clc;clear;
%数据输入
load runoff.txt;
load temperature.txt;
load precipitation.txt;
YearBegin=1957;
YearEnd=2008;
Years_pre=0;
k=3;
%调用函数过程
[X10,X10_pre,HCS1,MeanAE,MeanRE,AIC,NASH]=GM13Prediction(runoff,temperature,precipitation,YearBegin,YearEnd,Years_pre,k)