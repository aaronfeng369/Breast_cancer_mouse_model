load('C:\matlab\20160831\sample-2\relaxion10.mat')
addpath('C:\matlab')
t=[1:1:14000]';%总时间
c=[5 -1 -1 3 14000 14000 14000];%待求参数
%c=[4.908 -3.642 -0.994 3.538 3.851e3 39.204 3.851e3];
% P=[1,0.6009,0.1690,0.8031,10.3300,];%test
% data=P(1)-(P(2)*(1-exp((-1*t)/P(4)))+P(3)*(1-exp((-1*t)/P(5))));
data=Y;%所测数据
options=optimset('fminsearch');%调整拟合参数
options.TolX=1e-100;
options.Display='off';
options.MaxIter=10000000;
[x,sfval,sexit,soutput]=fminsearch(@fun,c,options,t,data);%拟合
for j=1:1:10%多次调整拟合结果
    c=x;
    [x,sfval,sexit,soutput]=fminsearch(@fun,c,options,t,data);
end
for i=1:1:14000
      G(i)=x(1)-(x(2)*(1-exp((-1*i)/x(5)))+x(3)*(1-exp((-1*i)/x(6)))+x(4)*(1-exp((-1*i)/x(7))));
end
G=G';
plot(t,data,'c',t,G,'k')
% subplot(2,1,1)
% plot(t,data,'b')
% subplot(2,1,2)
% plot(t,G,'r')
