function[beta_hat_SCC1,MSE_AVE,MSE_SIN,H]=ISCC_spatial_regressiongai_base(range,x,y,lon,lat,beta,sim_num,options)

%Estimate the regression coefficients using the SCC method 使用传统的SCC模型方法来模拟或实践 
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

% OUTPUT ARGUMENTS
% beta_hat_SCC1:  [sim_num,n,p] array storing the estimates of 
%                 regression coefficient
% MSE_AVE:        [sim_num,1] matrix storing the mean square error 
%                 of all explanatory variables
% MSE_SIN:        [sim_num,p] matrix storing the mean square error 
%                 for each explanatory variable
% H:              [np,np] matrix storing the location information of 
%                 the minimum spanning tree
[beta_hat_SCC1,MSE,H0]=SCC_spatial_regressiongai(x(range,:,:),y(range,:),lon,lat,beta,sim_num,[]);
MSE_SCC(:,1)=mean(MSE,2);%对应模拟次数下的三个系数的均方误差
MSE_AVE = mean(MSE,2);
MSE_SIN = MSE;
H = H0;
end