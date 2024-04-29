% Indentation fitting using modified prony series
% Ref: Qiang, B., Greenleaf, J., Oyen, M., & Zhang, X. (2011).  IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 58(7), 1418�C1429. http://doi.org/10.1109/TUFFC.2011.1961
%       Wang, X., Schoen, J. A., & Rentschler, M. E. (2013).  Journal of the Mechanical Behavior of Biomedical Materials, 20, 126�C136. http://doi.org/10.1016/j.jmbbm.2013.01.007
%
% Output Variables:
% aa - number of whole test data
% thk - thicknesses of samples
% meanthk - mean value of thk
% datalab0 - data acquired when in initial contact�սӴ�ʱ����
% datalab1 - actual indentation data ��ѹ��ȫ������ [180000points,force-displacement,5depth,6samples]
% indenterA��indenter areaѹͷ���
% indenterR��indenter radiusѹͷ�뾶
% [maxF, maxF_loc] - maximum value and index of the orignial force data acquired ԭ�������ֵ
% [minF, minF_loc] - minimum value and index of the orignial force data acquired ԭ������Сֵ
% [max2, max2loc] - maximum value and index of the flipped
%                   indentation-relaxation curve
% pts - [1x5], max2loc index value for all the indentation depth curves����ѹ�κͻָ��ηֽ��
% start_pt = 500;% the initial 1s (1e3 points) were acquired for waiting for the motor to warm up
% relaxFall - cell array containing all relaxation force (mN) [5depth,6samples]
% rampFall - cell array containing all ramp force (mN) [5depth,6samples]
% normFrelax - {5depth,6samples} normalized relaxation force
% normFramp - {5depth,6samples} normalized ramp force
% meanforce -  ��һ�������ѹ��ȵ�ƽ������
%
% Record of Revision
% Feb-17-2017===SQ,YF===Original Code
% May-17-2017===YF===Add comments and save individual f-d data

addpath('G:\Suda Research\Indentation device\mfiles\prony_v2')

%% clear and reset data ��������
clear all
clc
close
aa = 0;
addpath('C:\matlab\2017')
loadCellRatio = 9.8; % mN/V
laserSensorRatio = 6; % mm/V
indenterR = 1; % mm
indenterA = indenterR^2*pi; % mm2

%% read in data 
aa = aa+1;   % counting ����
fname0 = {'bc'};%�Ӵ�������
fname1 = {'bc_02','bc_04','bc_06','bc_08','bc_10'};%��ѹ������
color = {'r', 'y', 'g', 'c', 'b', 'm'};%��ɫ
color2 = {'m','r', 'y', 'g', 'c', 'b'};
% read in data
temp = sprintf([fname0{1}, '_%d.xlsx'], 0);%��ȡ�ļ���
temp1 = xlsread(temp);
windowSize = 100;
temp1 = filter(ones(1,windowSize)/windowSize,1,temp1); %��filter����ͼ��
temp1(:,[1,2])=temp1(:,[2,1]);
if length(temp1(:,1)) < 180000
    temp1(length(temp1(:,1))+1:180000,1)= temp1(end,1);
    temp1(length(temp1(:,1))+1:180000,2)= temp1(end,2);
end
datalab0(:,:,aa) = temp1;
inf = xlsread('breastca information.xlsx');% ��ȡ��ѹ���
thk(aa) = inf(5,1);
for ii =1:5;
    %  indentation data
    temp2 = 0;
    temp = sprintf([fname1{ii}, '_%d.xlsx'], 1);
    temp2 = xlsread(temp);
    temp2 = filter(ones(1,windowSize)/windowSize,1,temp2); %��filter����ͼ��
    temp2(:,[1,2])=temp2(:,[2,1]);
    if length(temp2(:,1)) < 180000
        temp2(length(temp2(:,1))+1:180000,1)= temp2(end,1);
        temp2(length(temp2(:,1))+1:180000,2)= temp2(end,2);
    end
    datalab1(:,:,ii,aa) = temp2;
end
%save('C:\\matlab\\2017\\Ca_mean\\MeanCa_patient.mat')%����ȫ������

%% grouping and normalization
meanthk = mean(thk);%����ƽ����ѹ���
start_pt = 500;% the initial 1s (1e3 points) were acquired for waiting for the motor to warm up
smoothRatio = 2000;

for jj = 1:size(datalab1, 3)
    for ii = 1:size(datalab1, 4)
        
        % get one FD curve
        data1 = datalab1(:,:,jj,ii);
        
        % indentation force section
        [maxD, maxD_loc] = max(data1(start_pt:end, 2));
        contactPt = start_pt + maxD_loc - 1;
         % find the minimal indentation force and its location, this is the peak force
        [minF, minF_loc] = min(data1(contactPt:end,1));
        % the minimal force indentation location index in the whole vector 
        peakPt = minF_loc + contactPt - 1;
        
        rampF = data1(contactPt:peakPt,1);
        rampDisp = data1((contactPt:peakPt),2); 
        relaxF = data1(peakPt:end, 1);
        relaxDisp = data1(peakPt:end, 2);
        contactPts(jj,ii) = contactPt;
        peakPts(jj,ii) = peakPt;
        
        % find the max indentation force/displacement and its location, this is the original data before flipping 
        [maxF, maxF_loc] = max(data1(contactPt:end, 1));
        
        % ramp force/displacement - corrected
        rampforce = (-rampF + maxF)*loadCellRatio;
        rampDisp_cor = (rampDisp(1) - rampDisp)*laserSensorRatio;
        % relaxation force/displacement - corrected
        relaxationforce = (-relaxF + maxF)*loadCellRatio; 
        relaxDisp_cor = (rampDisp(1) - relaxDisp)*laserSensorRatio;
        
        % smooth data
        relaxationforce = smooth(relaxationforce, smoothRatio);
        rampforce = smooth(rampforce, smoothRatio);
        relaxDisp_cor = smooth(relaxDisp_cor, smoothRatio);
        rampDisp_cor = smooth(rampDisp_cor, smoothRatio);
        
        % ramp and relaxation force after y axis correction
        relaxationforce = relaxationforce - rampforce(1);% ����0��ʼ
        rampforce = rampforce - rampforce(1);
        relaxFall{jj,ii} = relaxationforce;
        rampFall{jj,ii} = rampforce;
        Fall{jj,ii} = [rampforce; relaxationforce];
        
        % normalization with respect to Fmax
        normFrelax{jj,ii} = relaxationforce./relaxationforce(1); %���ݹ�һ��
        normFramp{jj,ii} = rampforce./relaxationforce(1);
        normFall{jj,ii} =  [rampforce; relaxationforce]./relaxationforce(1);
        
        % for displacement
        rampDispAll{jj,ii} = rampDisp_cor;
        relaxDispAll{jj,ii} = relaxDisp_cor;
        dispAll{jj,ii} = [rampDisp_cor; relaxDisp_cor];
    end
end

% check data
for ii = 1:size(relaxFall,2)
    figure; % for each indentation depth
    for jj = 1:size(relaxFall,1)
        subplot(1,2,1); plot(Fall{jj,ii},color{jj}); hold on
        subplot(1,2,2); plot(dispAll{jj,ii}, color{jj})
%         subplot(1,2,1); plot(rampDispAll{jj,ii},color{jj}); hold on
%         subplot(1,2,2); plot(relaxDispAll{jj,ii}, color{jj})
    end
end

% legend('2%','4%','6%','8%','10%','location','NorthEast');%���ע��
% ylabel('F(mN)','FontName','Times New Roman','FontSize',24);%����y��
% xlabel('T(ms)','FontName','Times New Roman','FontSize',24);%����x��
% set(gca,'FontName','Times New Roman','FontSize',24,'LineWidth',1.5);%�����������߿�
% title('F-T','FontName','Times New Roman','FontSize',28);%��ӱ���
%
% fh2 = figure;
% plot(max2loc(:), max2(:), 'o-k', 'LineWidth',4);%������ߵ������ͼ
% ylabel('F(mN)','FontName','Times New Roman','FontSize',24);%����y��
% xlabel('T(ms)','FontName','Times New Roman','FontSize',24);%����x��
% set(gca,'FontName','Times New Roman','FontSize',24,'LineWidth',1.5);%�����������߿�
% title('top point line','FontName','Times New Roman','FontSize',28);%��ӱ���

% axis square;
% print(fh1, '-dtiff', '-r300','F-T');%����ͼƬ
% print(fh2, '-dtiff', '-r300','top point line');

%% plot meandata f-t figure ��ƽ��f-tͼ��
start_pt = 500;% the initial 1s (1e3 points) were acquired for waiting for the motor to warm up
for j = 1:5 %�����ѹ��ȵ�ƽ��ֵ
    data1 = 0;
    data1 = normFall(:,:,j,1);
    meanforce(j,:) = data1;%��ȡ����һ������
    for jj = 2:aa%֮���������������������
        temp3 = normFall(:,:,j,jj);
        meanforce(j,:) = meanforce(j,:) + temp3;
    end
end
meanforce = meanforce./aa;%���ƽ����
fh1 = figure;
for ii=1:5 %��ѹ����ָ��ηֿ�
    data1 = 0;
    data1 = meanforce(ii,:);
    meantlen(ii) = length(data1);% ���ݳ���
    
    % find the max indentation force and its location, this is the original data before flipping ���ֵ������֮��ͼ��ת
    [meanmaxF(ii) ,meanmaxF_loc(ii)] = max(data1);
    % find the minimal indentation force and its location ��Сֵ
    [meanminF(ii), meanminF_loc(ii)] = min(data1);
    % the minial force indentation location index in the whole vector ��Сֵ��λ�ã�Ҳ����ѹ�ͻָ��ķֽ��
    peakPt(ii) = meanmaxF_loc(ii) - 1;
    % relaxation force section �ָ���������
    meanrelaxF = data1((peakPt(ii):meantlen(ii)));
    % indentation force section ��ѹ��������
    meanrampF = data1((1:peakPt(ii)));
    
    %     relaxationforce = smooth(relaxationforce,2000);
    %     rampforce = smooth(rampforce,2000);
    
    % ramp and relaxation force after y axis correction
    meanrelaxForce{ii} =  meanrelaxF;
    meanrampForce{ii} = meanrampF;
    
    % vector length of the ramp and relaxation sections
    meanreltimelen(ii) = length(meanrelaxF);%��ȡ�������ݳ���
    meanramplen(ii) = length(meanrampF);
    % assmble two sections
    %     plot(meanforce,color{ii});%��ͼ
    %     hold on
    plot((1/1000:0.001:meanreltimelen(ii)/1000),meanrelaxF,color{ii});%��ͼ
    hold on
    %         plot(meanramplen(ii):(meanreltimelen(ii)+meanramplen(ii)-1),(relaxationforce),color2{ii});
    %         hold on
    %data1(:,:) = 0;
end
box off
h1 = legend('2%','4%','6%','8%','10%','location','NorthEast');%���ע��
set(h1,'Box','off');
axis square;
ylabel('normalized stress','FontName','Times New Roman','FontSize',24);%����y��
xlabel('T (s)','FontName','Times New Roman','FontSize',24);%����x��
set(gca,'FontName','Times New Roman','FontSize',24,'LineWidth',1.5);%�����������߿�
%title('F-T','FontName','Times New Roman','FontSize',28);%��ӱ���
% ������ߵ�����
% fh2 = figure;
% plot(max2loc(:), max2(:), 'o-k', 'LineWidth',4);%������ߵ������ͼ
% ylabel('F(mN)','FontName','Times New Roman','FontSize',24);%����y��
% xlabel('T(ms)','FontName','Times New Roman','FontSize',24);%����x��
% set(gca,'FontName','Times New Roman','FontSize',24,'LineWidth',1.5);%�����������߿�
% title('top point line','FontName','Times New Roman','FontSize',28);%��ӱ���
save meanCancerData_patient.mat%����ȫ������

%% ��������ͼ
fh2 = figure;
for jj = 5:5%�����������ݸ��Ե�ƽ��ֵ����Ѱ�����ֵ��Сֵ
    for ii = 1:aa
        comparedata(jj,ii) = mean(normFrelax{jj,ii});
        %         plot((1/1000:0.001:length(normFrelax{jj,ii})/1000),normFrelax{jj,ii});
        %         hold on
    end
    hold on
    %����ƽ��������
    plot((1/1000:0.001:meanreltimelen(jj)/1000),meanrelaxForce{jj},'k','LineWidth',3.5)
    hold on
    %�ҵ������С����
    [mindata(jj),min_loc(jj)] = min(comparedata(jj,1:end));
    [maxdata(jj),max_loc(jj)] = max(comparedata(jj,1:end));
    minlength = length(normFrelax{jj,min_loc(jj)});
    maxlength = length(normFrelax{jj,max_loc(jj)});
    plot((1/1000:0.001:minlength/1000),normFrelax{jj,min_loc(jj)},'k','LineWidth',1)%��ͼ
    hold on
    plot((1/1000:0.001:maxlength/1000),normFrelax{jj,max_loc(jj)},'k','LineWidth',1)
    hold on
    %ʱ��msת��Ϊs
    x_minrange = 1/1000:0.001:minlength/1000;
    x_maxrange = 1/1000:0.001:maxlength/1000;
    y1 = normFrelax{jj,max_loc(jj)};  % ���������ϱ߽�
    y2 = normFrelax{jj,min_loc(jj)};  % ���������±߽�
    y_upper = y1(1:maxlength)';
    y_lower = y2(1:minlength)';
    % �����м�����
    patch([x_minrange, fliplr(x_maxrange)], [y_lower, fliplr(y_upper)], 'r')
    alpha(0.5) % ����͸����
end
box off
grid on
%axis square;
%legend('mean','min','max','location','NorthEast');%���ע��
ylabel('normalized stress','FontName','Times New Roman','FontSize',24);%����y��
xlabel('T (s)','FontName','Times New Roman','FontSize',24);%����x��
set(gca,'FontName','Times New Roman','FontSize',24,'LineWidth',1.5);%�����������߿�
%title('4%','FontName','Times New Roman','FontSize',28);%��ӱ���
print(fh2, '-dtiff', '-r300','patient_10%');%����ͼ��

