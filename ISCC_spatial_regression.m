function [beta_hat_SCC1_end,MSE_SCC,MSE_MEAN,A_SCC]=ISCC_spatial_regression(range,x,y,...
                                           lon,lat,beta,sim_num,options,zhibiao,c,p,diedai)
%Estimate the regression coefficients using the ISCC method
%使用ISCC模型方法来模拟或实践(固定p为3)

% Input ARGUMENTS
% sim_num: number of simulations 
% range:   The range of simulations in x and y (ex: range = 1:sim_num)
% x:       [sim_num,n,p] array storing the p explanatory variables on n locations
% y:       [sim_num,n] matrix storing the response variables on n locations
% lon:     [n,1] matrix storing the longitudes of n locations
% lat:     [n,1] matrix storing the latitudes of n locations
% beta:    [n,p] matrix storing the true values of regression coefficients
% options: a structure setting the parameters used in the SCC
%          see function SCC.m for details
% zhibiao: threshold for weight adjustment
% c:       factor of zhibiao change
% diedai:  number of steps of ISCC model iteration
% p:       number of variables

% OUTPUT ARGUMENTS
% beta_hat_SCC1_end:  [sim_num,n,p] array storing the estimates of 
%                      regression coefficient
% MSE_SCC:        [sim_num,diedai+1] matrix storing the mean square error 
%                 of each simulation in each iteration of the ISCC model
% MSE_MEAN:       [1,diedai+1] matrix storing the mean square error 
%                 for each explanatory variable in each iteration of the ISCC model
% A_SCC:          [sim_num,diedai+1] matrix storing the A value of each simulation 
%                 in each iteration of the ISCC model
                 

[beta_hat_SCC1,MSE_AVE,~,H]=ISCC_spatial_regression_base...
                           (range,x,y,lon,lat,beta,sim_num,options);

H0=H;
Hqian=H0; 
dzengjiahe1=[];
dzengjiahe2=[];
dzengjiahe3=[];
beta_hat_SCC1_last=beta_hat_SCC1;%初始估计参数传参过渡 
beta_hat_SCC1_ori=beta_hat_SCC1;%初始估计参数传参过渡 (用这个画图)
beta_hat_SCC1_end=zeros(sim_num,1000,3);%存储估计出的参数
A_SCC=nan(sim_num,diedai+1);%存储真实边界的边数
MSE_SCC=nan(sim_num,diedai+1);
MSE_SCC(:,1)=MSE_AVE;

for j=1:sim_num %一个一个做范围指标
    Hqian=H0; 
    dzengjiahe1=[];
    dzengjiahe2=[];
    dzengjiahe3=[];
    beta_hat_SCC1=beta_hat_SCC1_ori(j,:,:);
    beta_hat_SCC1_last=beta_hat_SCC1_ori(j,:,:);
    change = zeros((diedai-1),1);%迭代两步间的系数变动
for i=1:diedai
    if i>=2
    change(i-1,1) = sum(sum(abs(beta_hat_SCC1-beta_hat_SCC1_last)));
    end
%    j可显示模拟进程
%    i可显示模拟进程
    beta_hat_SCC1_last = beta_hat_SCC1;%用来计算两步间系数变动差的中间变量
    [A_SCC(j,1),A_SCC(j,i+1)]=A(lon,lat,beta,zhibiao,p,0.1,Hqian,...
                              beta_hat_SCC1,dzengjiahe1,dzengjiahe2,dzengjiahe3);
    if i>=2
    panduan_change = change(i-1,1);
    if panduan_change==0
        zhibiao=1e-3;
    else
        zhibiao=zhibiao*c;
    end
    end
    [beta_hat_SCC2,MSE,Hhou,dzengjiahe1hou,dzengjiahe2hou,dzengjiahe3hou]=...
                                SCC_spatial_regressiongai2fenpro(x(range(j),:,:),y(range(j),:),...
                                lon,lat,beta,1,[],beta_hat_SCC1,diedai,Hqian,dzengjiahe1,...
                                dzengjiahe2,dzengjiahe3,zhibiao);
    MSE_SCC(j,i+1)=mean(MSE,2);
    beta_hat_SCC1=beta_hat_SCC2;
    Hqian=Hhou;
    dzengjiahe1=dzengjiahe1hou;
    dzengjiahe2=dzengjiahe2hou;
    dzengjiahe3=dzengjiahe3hou;
end
beta_hat_SCC1_end(j,:,:)=beta_hat_SCC1;%保存估计出的参数
end
hengzhou=0:diedai;
MSE_MEAN=mean(MSE_SCC);
end