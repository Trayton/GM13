%-----------------------------------------------------------------
%GM(1,3)ģ��
%��ֱ������GM13_Main.m���ɡ�
%ԭ�������ۼ�����X1��ͨ����С���˷����ƻ�΢�ַ��̵Ĳ���������������׻���̬���̻�ԭX0��
%�����������ѧ������1990�����ĵ˾�������ά��ɫ�滮����һ��GM��1��N)ģ�ͣ�
%��������������
%runoff.txtΪ��Ԥ������runoff�������е�ʱ���������ݣ�
%temperature.txtΪӰ������temperature�������е�ʱ���������ݣ�
%precipitation.txtΪӰ������precipitation�������е�ʱ���������ݣ�
% YearBeginΪ��ʼ��ݣ�
% YearEndΪ��ֹ��ݣ�
% Years_preΪ��Ԥ���������
% kΪGM13�Ĳ�����������Ȼ��ͬ��GM13�е�3��
%�������������
% X10Ϊ��Ԥ������runoff��ʵ�ʹ۲�ֵ��
% X10_preΪ��Ԥ������runoff��ģ��ֵ��
% HCS1ΪGM13ģ�͵Ĳ���a��b1��b2���ֱ���runoff��temperature��precipitation�Ĳ�����
% MeanAE����˼����ƽ���������
% MeanRE��ƽ��������
% AICΪ��С��Ϣ׼���ǲ�������k��Ԥ��в�ĺ�����ԽСԽ�ã�
% NASHϵ��Ϊˮ��ѧ��������ģ��ȷ���Ե�ϵ����
% Author: Chong Wang E-mail:chong0314@gmail.com
%------------------------------------------------------------------
clc;clear;
%��������
load runoff.txt;
load temperature.txt;
load precipitation.txt;
YearBegin=1957;
YearEnd=2008;
Years_pre=0;
k=3;
%���ú�������
[X10,X10_pre,HCS1,MeanAE,MeanRE,AIC,NASH]=GM13Prediction(runoff,temperature,precipitation,YearBegin,YearEnd,Years_pre,k)