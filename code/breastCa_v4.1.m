% Indentation fitting using modified prony series
% Ref: Qiang, B., Greenleaf, J., Oyen, M., & Zhang, X. (2011).  IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 58(7), 1418–1429. http://doi.org/10.1109/TUFFC.2011.1961
%       Wang, X., Schoen, J. A., & Rentschler, M. E. (2013).  Journal of the Mechanical Behavior of Biomedical Materials, 20, 126–136. http://doi.org/10.1016/j.jmbbm.2013.01.007
% 
% Output Variables:
% RMSE - Root mean square error, objective function for optimization：拟合数据误差
% coe - [5x7], 拟合系数
% thk - thickness of the sample 样本厚度
% data0 - data acquired when in initial contact,[data points; force,displacment; 2,4,6,8,10% strain] 刚接触时数据
% data1 - actual indentation data, [data points; force,displacment; 2,4,6,8,10% strain] 下压后全部数据
% data2 - indentation data with force filtered by an moving average filter
% indenterA - indenter area (mm2)压头面积
% indenterR - indenter radius (mm)压头半径
% k - slope of the linear fitting of the indentation ramp 压下段线性拟合斜率
% [maxF, maxF_loc] - maximum value and index of the orignial force data acquired 原数据最大值
% [minF, minF_loc] - minimum value and index of the orignial force data acquired 原数据最小值
% [max2, max2loc] - maximum value and index of the flipped
%                   indentation-relaxation curve
% pts - [1x5], max2loc index value for all the indentation depth curves：下压段和恢复段分界点
% ramplen - vector length of the indentation ramp. Time in ms下压段数据长度
% reltimelen - vector length of the relaxation section. For fitting, it represents time in ms. 恢复段数据长度
% relaxForce{ii} - cell array containing all relaxation force;
% rampForce{ii} - cell array containing all ramp force;
% totalF{ii} - cell array containg all total force
% totalDisp - all displacement data
% Units: kPa, s*Pa,
% 
% Record of Revision
% Feb-17-2017===SQ,YF===Original Code
% Mar-04-2017===YF===Modify fitting to use RMSE, add comments
% Mar-27-2017===YF===Organize force and displacement data
% Mar-29-2017===YF===smooth data before analysis

addpath('G:\Suda Research\Indentation device\mfiles\prony_v2')

%% read in data 读取数据
clear all
clc
close 
addpath('C:\matlab\2017')
fname0 = {'bc'};%接触数据名
fname1 = {'bc_02','bc_04','bc_06','bc_08','bc_10'};%下压数据名
color = {'r', 'y', 'g', 'c', 'b', 'm'};%颜色
loadCellRatio = 27; % mN/V
laserSensorRatio = 6; % mm/V
indenterW = 5.79*9.8; % mg*g=mN, this is the indenter head weight before contact
indenterR = 1; % mm
indenterA = indenterR^2*pi; % mm2
%load('thick.mat'); % sample thickness 样本厚度
%thk = thickness;
thk = 4.63;% sample thickness
% read in data
temp = sprintf([fname0{1}, '_%d.mat'], 0);%读取文件名
temp1 = load(temp);
data0(:,:) = temp1.Y;
for ii =1:5;
    %  indentation data
    temp = sprintf([fname1{ii}, '_%d.mat'], 1);
    temp1 = load(temp);
    data1(:,:, ii) = temp1.Y;
end
   
%% plot f-t figure, find the rule of top points 画全部f-t图像并比较最高点的关系
start_pt = 1000;% the initial 1s (1e3 points) were acquired for waiting for the motor to warm up
for ii = 1:5
    data2(:,1,ii) = smooth(data1(:, 1, ii), 300);
    data2(:,2,ii) = smooth(data1(:,2,ii), 300);
end

fh1 = figure;
for ii=1:5
    % temp = 
    tlen(ii) = length(data2(:, 1, ii));% 数据长度
    
    % find the minimal indentation force and its location 最小值，前1000个点无用可删除
    [minF(ii), minF_loc(ii)] = min(data2(start_pt:tlen(ii), 1, ii));
    % the minial force indentation location index in the whole vector 最小值的位置，也是下压和恢复的分界点
    pts(ii) = minF_loc(ii) + start_pt - 1; 
     % relaxation force section 恢复部分数据
    relaxF = data2((pts(ii):tlen(ii)), 1, ii);
    % indentation force section 下压部分数据
    rampF = data2((start_pt:pts(ii)),1,ii);
    
    % find the max indentation force and its location, this is the original data before flipping 最大值，便于之后图像翻转
    [maxF(ii) ,maxF_loc(ii)] = max(data2(start_pt:tlen(ii), 1, ii));
    % relaxation force - corrected
    relaxationforce = (((relaxF-minF(ii))*(-1)) + (maxF(ii)-minF(ii)))*loadCellRatio; % 转换为力，并将图像上下翻转
    % ramp force - corrected
    rampforce = (((rampF - minF(ii))*(-1))+(maxF(ii)-minF(ii)))*loadCellRatio;
%     % smooth the data and reduce noise
%     relaxationforce = smooth(relaxationforce,2000);
%     rampforce = smooth(rampforce,2000);
    % ramp and relaxation force after y axis correction
    relaxationforce = relaxationforce - rampforce(1);% 起点从0开始
    rampforce = rampforce - rampforce(1);
    relaxForce{ii} = relaxationforce;
    rampForce{ii} = rampforce;
    
    % vector length of the ramp and relaxation sections
    reltimelen(ii) = length(relaxF);%读取各端数据长度
    ramplen(ii) = length(rampF);
        % overall curve max point
    max2(ii) = relaxationforce(1);%读取各数据最高点数值和位置
    max2loc(ii) = ramplen(ii);
    
    % assmble two sections
    force(1:ramplen(ii)) = rampforce;%两端数据合并，便于注释
    force((ramplen(ii)+1):(reltimelen(ii)+ramplen(ii))) = relaxationforce;
    totalF{ii} = force;
    
    plot(force,color{ii});%画图
    hold on
    %     plot((rampforce),color{ii});
    %     hold on
    %     plot(ramplen(ii):(reltimelen(ii)+ramplen(ii)-1),(relaxationforce),color{ii});
    %     hold on
end

legend('2%','4%','6%','8%','10%','location','NorthEast');%添加注释
ylabel('F(mN)','FontName','Times New Roman','FontSize',24);%设置y轴
xlabel('T(ms)','FontName','Times New Roman','FontSize',24);%设置x轴
set(gca,'FontName','Times New Roman','FontSize',24,'LineWidth',1.5);%设置坐标轴线宽
title('F-T','FontName','Times New Roman','FontSize',28);%添加标题

fh2 = figure;
plot(max2loc(:), max2(:), 'o-k', 'LineWidth',4);%画出最高点的连线图
ylabel('F(mN)','FontName','Times New Roman','FontSize',24);%设置y轴
xlabel('T(ms)','FontName','Times New Roman','FontSize',24);%设置x轴
set(gca,'FontName','Times New Roman','FontSize',24,'LineWidth',1.5);%设置坐标轴线宽
title('top point line','FontName','Times New Roman','FontSize',28);%添加标题

% axis square;
% print(fh1, '-dtiff', '-r300','F-T');%保存图片
% print(fh2, '-dtiff', '-r300','top point line');

%% find f-d figures and k values 画f-d图像并计算各端的k值
% save all displacement data
totalDisp = squeeze(data2(:,2,:));

fh2 = figure;
del_pts = 30; % time difference between force and displacement 位移与压力传感器数据的时间差
for ii = 2:5
    y = polyfit(((data2((start_pt+del_pts):(ramplen(ii)+start_pt+del_pts-1), 2, ii) - data2((start_pt+del_pts),2, ii))*(-1*laserSensorRatio))...
        ,(((data2(start_pt:(ramplen(ii)+start_pt-1),1, ii) - data2(start_pt,1, ii))*(-1*loadCellRatio))),1);%对下压端的翻转图像进行直线拟合
    k(ii) = y(1);%保存k值
    yy = polyval(y,((data2((start_pt+del_pts):ramplen(ii)+(ramplen(ii)+start_pt+del_pts-1),2, ii)...
        -data2((start_pt+del_pts),2, ii))*(-1*laserSensorRatio)));%求出拟合线
    displacement(ii) = max((data2((start_pt+del_pts):(ramplen(ii)+start_pt+del_pts-1),2, ii)...
        -data2((start_pt+del_pts),2, ii))*(-1*laserSensorRatio));%位移距离，mm
    plot(((data2((start_pt+del_pts):(ramplen(ii)+start_pt+del_pts-1),2, ii)-data2((start_pt+del_pts),2, ii))*(-1*laserSensorRatio)), ...
        (((data2(start_pt:(ramplen(ii)+start_pt-1),1, ii)-data2(start_pt,1, ii))*(-1*loadCellRatio))), color{ii}); %画图
    hold on
    %     plot(((data2(1000:ramplen(ii)+999,2, ii)-data2(1000,2, ii))*(-1*laserSensorRatio)),yy,color{ii});
    %     hold on
end
legend('4%','6%','8%','10%','location','SouthEast');%设置图像格式
ylabel('F(mN)','FontName','Times New Roman','FontSize',24);
xlabel('d(mm)','FontName','Times New Roman','FontSize',24);
set(gca,'FontName','Times New Roman','FontSize',24,'LineWidth',1.5);
title('F-d','FontName','Times New Roman','FontSize',28);
% print(fh2, '-dtiff', '-r300','F-D');%保存图像

fh3 = figure;
for ii = 2:5 %画出2s的位移实际数据
    plot((data2((start_pt+del_pts):ramplen(ii)+1999,2, ii)-data2((start_pt+del_pts),2, ii))*(-1*laserSensorRatio),color{ii});
    hold on
end
hold on
for ii = 2:5 %画出计算位移距离所用数据长度，比较与实际数据长度之间的差异
    plot((data2((start_pt+del_pts):(ramplen(ii)+start_pt+del_pts-1),2, ii)-data2((start_pt+del_pts),2, ii))*(-1*laserSensorRatio),'k');
    hold on
end
legend('4%','6%','8%','10%','location','SouthEast');%设置图像格式
ylabel('d(mm)','FontName','Times New Roman','FontSize',24);
xlabel('T(ms)','FontName','Times New Roman','FontSize',24);
set(gca,'FontName','Times New Roman','FontSize',24,'LineWidth',1.5);
title('d-T','FontName','Times New Roman','FontSize',28);
%save CancerData.mat

%% fitting relaxation curve to prony series

termOpt = 3; % term of prony series 阶数
N = 3; % number of cycles of circuit fitting

fh3 = figure;
for ii = 2:5
    %c1(ii,:)=[6 4 5 4 80];
    %c1(ii,:)=[8 5 3.5 34000 690000];
    %拟合初值
    %     c1(ii,:) = [13 5 5300]; % 1 term
    %     c2(ii,:) = [13 5 5 2300 33000]; % 2 term
    %     %c2(ii,:) = [13 5 5.2 5300 120000]; % 2 term 1min
    %     %c3(ii,:) = [13 5 55 58 2300 33000 160000];  % 3 term 3min 10%
    %     %c3(ii,:) = [10 5 10 35 2300 33000 160000]; % 3 term 3min 8% 5mm mouse
    %     c3(2,:) = [0.4 0.3 0.4 0.4 2300 33000 160000]; % 3 term 3min 4% 2mm
    %     c3(3,:) = [0.05 0.1 0.05 0.1 2300 33000 160000]; % 3 term 3min 6% 2mm
    %     %c3(4,:) = [0.1 0.3 0.1 0.3 2300 33000 160000]; % 3 term 3min 8% 2mm
    %     c3(4,:) = [0.1 0.1 0.1 0.1 300 3000 16000];
    %     c3(5,:) = [0.1 0.1 0.1 0.1 2300 33000 160000]; % 3 term 3min 10% 2mm
    %
    %     c3(2,:) = [0.4 0.3 0.4 0.6 2300 33000 160000]; % 3 term 3min 4% 5mm
    %     c3(3,:) = [0.5 0.3 0.4 0.5 2300 33000 160000]; % 3 term 3min 6% 5mm
    %     c3(4,:) = [0.5 0.3 0.4 0.8 2300 33000 160000]; % 3 term 3min 8% 5mm
    %     c3(5,:) = [0.5 0.3 0.4 1 2300 33000 160000]; % 3 term 3min 10% 5mm
    %c3(ii,:) = [23 5 65 68 2300 33000 160000];  % 3 term 3min 8% 5mm
    %c3(ii,:) = [3 5 55.2 5 900 9000 140000]; % 3 term 1min
    %c3(ii,:) = [12.7 3 3.4 5 8000 960 140000];
    c = 0;
    x = 0;
    t1 = [1:1:reltimelen(ii)]';%恢复时间
    options1 = optimset('fminsearch');%拟合选项设置
    options1.TolX = 1e-100;
    options1.TolFun = 1e-100;
    options1.MaxIter = 1000000000;
    options1.Display = 'off';
    
    options = optimset('fmincon');%拟合选项设置
    options.TolX = 1e-100;
    options.TolFun = 1e-100;
    options.MaxIter = 1000000000;
    options.Display = 'off';
    
    options2=optimset('ga');%调整拟合参数
    options2.TolX=1e-100;
    options2.TolCon=1e-100;
    options2.TolFun=1e-100;
    options2.MaxIter=1000000000;
    options2.Display='off';
    
    relaxationforce =  smooth(relaxForce{ii},500);%数据读取，加数据降噪处理
    %     if termOpt == 1
    %         % 1 term
    %         c = c1(ii,:);
    %     elseif termOpt == 2
    %         % 2 term
    %         c = c2(ii,:);
    %     else
    %         % 3 term
    %         c = c3(ii,:);
    %     end
    
    %norelforce=relaxationforce./max2(ii);%归一化
    p = sqrt(indenterR*displacement(ii))/thk;
    omeg = (1+1.33*p+1.283*p^2+0.769*p^3+0.0975*p^4);
    rate = (8*indenterR*displacement(ii)*omeg)/3;%求常数项
    %[x,sfval1,sexit1,soutput1] = fminsearch(@pronySeries, c, options1, t1, rate, relaxationforce,termOpt);%拟合
    [x,sfval1,sexit1,soutput1]=ga(@(x)  pronySeries(x,t1,rate,relaxationforce,termOpt),7,[],[],[],[],[0 0 0 0 0 0 0],[],[],options2);%遗传算法拟合求大概值
    % circuit fitting in order to reduce errors
    for i = 1:N
        c = x;
        %[x,sfval1,sexit1,soutput1] = fminsearch(@pronySeries, c, options, t1, rate, relaxationforce,termOpt);%拟合
        [x,sfval1,sexit1,soutput1]=fmincon(@(x) pronySeries(x,t1,rate,relaxationforce,termOpt),c,[],[],[],[],[0.001 0.001 0.001 0.001 0 0 0],[],[],options);%多次拟合求出精确值
    end
    coe(ii,:) = x; %保存拟合系数
    Fcal = 0; %清零,防止数据长度不一
    
    % fitted force
    %求出拟合后的结果
    if termOpt == 1
        % 1 term
        for t=1:1:reltimelen(ii)
            Fcal(t) = rate*(x(1)+x(2)*exp(-t/(x(3)/x(2))));
        end
    elseif termOpt == 2
        % 2 terms
        for t=1:1:reltimelen(ii)
            Fcal(t) = rate*(x(1)+x(2)*exp(-t/(x(4)/x(2)))+x(3)*exp(-t/(x(5)/x(3))));
        end
    else
        % 3 terms
        for t=1:1:reltimelen(ii)
            Fcal(t) = rate*(x(1)+x(2)*exp(-t/(x(5)/x(2)))+x(3)*exp(-t/(x(6)/x(3)))+x(4)*exp(-t/(x(7)/x(4))));
        end
    end
    
    Fcal=Fcal';
    
    %MSE(ii) =  sum((relaxationforce-Fcal).^2)/reltimelen(ii);%求出误差系数
    RMSE(ii) = sqrt(sum((relaxationforce-Fcal).^2)/reltimelen(ii));
    xcoe(ii,1:7) = coe(ii,1:7); %合并数据，方便整理
    xcoe(ii,8) = RMSE(ii);
    
    %subplot(2,2,(ii-1))%画图
    plot(t1,relaxationforce,color{ii},'LineWidth',2);
    hold on
    plot (t1,Fcal,'k','LineWidth',2);
    
    %title(2*ii);
    
end
print(fh3, '-dtiff', '-r300','Fitting');%保存图像
save CancerData.mat%保存全部数据

% save individual data
fh5 = figure;
for ii = 1:5; plot(totalF{ii},color{ii}); hold on; end

for ii = 1:5; plot(data2(:,2,ii),color{ii}); hold on; end
force2 = totalF{1}; force4 = totalF{2}; force6 = totalF{3};  force8 = totalF{4};
save('indentF.mat', 'force2', 'force4', 'force6', 'force8')