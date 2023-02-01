function [count_yuanshi_z,count_hou_z]=A(lon,lat,beta,zhibiao,p,dc,Hqian,beta_hat_SCC1,dzengjiahe1,dzengjiahe2,dzengjiahe3)
% Calculate the A value of each iteration of the ISCC model
% 计算ISCC模型每步迭代的A值的大小（MST中跨越cluster的边数） 

% Input ARGUMENTS
% lon:             [n,1] matrix storing the longitudes of n locations
% lat:             [n,1] matrix storing the latitudes of n locations
% beta:            [n,p] matrix storing the true values of regression coefficients
% zhibiao          threshold for weight adjustment
% p:               number of variables
% dc:              a critical distance larger than which the two vertices will
%                  not be connected. This value is used only for the reduction 
%                  of computational burden.
% Hqian:           [np,np] matrix storing the location information of 
%                  the minimum spanning tree
% beta_hat_SCC1:   [1,n,p] array storing the estimates of 
%                  regression coefficient
% dzengjiahe1:     the position information of the previously broken edge 
%                  for the first explanatory variable
% dzengjiahe2:     The position information of the previously broken edge 
%                  for the second explanatory variable
% dzengjiahe3:     The position information of the previously broken edge 
%                  for the third explanatory variable

% OUTPUT ARGUMENTS
% count_yuanshi_z: calculate the number of edges across cluster in the 
%                  original minimum spanning tree
% count_hou_z:     calculate the number of edges across cluster in the 
%                  current minimum spanning tree


%需要输入load传入lon、lat、beta，需要zhibiao、p、Hqian、beta_hat_SCC1_end(1,:,:)
%load(location);%加载的数据地址
%zhibiao=36.1189;
%p=3;
%dc=0.1;

n=length(lon);
for j=1:p
if j==1
dbeta=squeeze(beta_hat_SCC1(1,:,:));
dbeta(:,j)=dbeta(:,j)/std(dbeta(:,j));
d=abs((dbeta(:,j)-dbeta(:,j)'));
[m,n]=size(d);  
Hqianbianliang=Hqian((n-1)*(j-1)+1:(n-1)*j,n*(j-1)+1:n*j);
jilu=[];
   for l=1:n-1
       jilu=[jilu;find(Hqianbianliang(l,:)~=0)];
   end
jilu;
dzengjia=[];
  if j==1
     for i=1:length(dzengjiahe1)
         d(dzengjiahe1(i,1),dzengjiahe1(i,2))=inf;
         d(dzengjiahe1(i,2),dzengjiahe1(i,1))=inf;
     end
  elseif j==2
     for i=1:length(dzengjiahe2)
         d(dzengjiahe2(i,1),dzengjiahe2(i,2))=inf;
         d(dzengjiahe2(i,2),dzengjiahe2(i,1))=inf;
     end
  elseif j==3
     for i=1:length(dzengjiahe3)
         d(dzengjiahe3(i,1),dzengjiahe3(i,2))=inf;
         d(dzengjiahe3(i,2),dzengjiahe3(i,1))=inf;
     end
  end
  
   for h=1:length(jilu)
       if abs(d(jilu(h,1),jilu(h,2)))<zhibiao  
           d(jilu(h,1),jilu(h,2))=1;
           d(jilu(h,2),jilu(h,1))=1;
       else
           d(jilu(h,1),jilu(h,2))=inf;  
           d(jilu(h,2),jilu(h,1))=inf;
           dzengjia=[dzengjia;jilu(h,1),jilu(h,2);jilu(h,2),jilu(h,1)];   
       end
   end
   if j==1
      dzengjiahe1=[dzengjiahe1;dzengjia];
   elseif j==2
      dzengjiahe2=[dzengjiahe2;dzengjia];
   elseif j==3
       dzengjiahe3=[dzengjiahe3;dzengjia];
   end
   d(find(d~=1&d~=inf))=1;
   dhalf=sqrt((lon-lon').^2+(lat-lat').^2);
   d1=d.*dhalf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if j==2 
dbeta=squeeze(beta_hat_SCC1(1,:,:));
dbeta(:,j)=dbeta(:,j)/std(dbeta(:,j));
d=abs((dbeta(:,j)-dbeta(:,j)'));
range=max(max(d))-min(min(d));
shangjie=999999;
[m,n]=size(d);  
Hqianbianliang=Hqian((n-1)*(j-1)+1:(n-1)*j,n*(j-1)+1:n*j);
jilu=[];
   for l=1:n-1
       jilu=[jilu;find(Hqianbianliang(l,:)~=0)];
   end
jilu;
dzengjia=[];
  if j==1
     for i=1:length(dzengjiahe1)
         d(dzengjiahe1(i,1),dzengjiahe1(i,2))=inf;
         d(dzengjiahe1(i,2),dzengjiahe1(i,1))=inf;
     end
  elseif j==2
     for i=1:length(dzengjiahe2)
         d(dzengjiahe2(i,1),dzengjiahe2(i,2))=inf;
         d(dzengjiahe2(i,2),dzengjiahe2(i,1))=inf;
     end
  elseif j==3
     for i=1:length(dzengjiahe3)
         d(dzengjiahe3(i,1),dzengjiahe3(i,2))=inf;
         d(dzengjiahe3(i,2),dzengjiahe3(i,1))=inf;
     end
  end
  
   for h=1:length(jilu)
       if abs(d(jilu(h,1),jilu(h,2)))<zhibiao  
           d(jilu(h,1),jilu(h,2))=1;
           d(jilu(h,2),jilu(h,1))=1;
       else
           d(jilu(h,1),jilu(h,2))=inf; 
           d(jilu(h,2),jilu(h,1))=inf;
           dzengjia=[dzengjia;jilu(h,1),jilu(h,2);jilu(h,2),jilu(h,1)];  
       end
   end
   if j==1
      dzengjiahe1=[dzengjiahe1;dzengjia];
   elseif j==2
      dzengjiahe2=[dzengjiahe2;dzengjia];
   elseif j==3
       dzengjiahe3=[dzengjiahe3;dzengjia];
   end
   d(find(d~=1&d~=inf))=1;
dhalf=sqrt((lon-lon').^2+(lat-lat').^2);
d2=d.*dhalf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if j==3
dbeta=squeeze(beta_hat_SCC1(1,:,:));
dbeta(:,j)=dbeta(:,j)/std(dbeta(:,j));
d=abs((dbeta(:,j)-dbeta(:,j)'));
range=max(max(d))-min(min(d));
shangjie=999999;
[m,n]=size(d);  
Hqianbianliang=Hqian((n-1)*(j-1)+1:(n-1)*j,n*(j-1)+1:n*j);
jilu=[];
   for l=1:n-1
       jilu=[jilu;find(Hqianbianliang(l,:)~=0)];
   end
jilu;
dzengjia=[];
  if j==1
     for i=1:length(dzengjiahe1)
         d(dzengjiahe1(i,1),dzengjiahe1(i,2))=inf;
         d(dzengjiahe1(i,2),dzengjiahe1(i,1))=inf;
     end
  elseif j==2
     for i=1:length(dzengjiahe2)
         d(dzengjiahe2(i,1),dzengjiahe2(i,2))=inf;
         d(dzengjiahe2(i,2),dzengjiahe2(i,1))=inf;
     end
  elseif j==3
     for i=1:length(dzengjiahe3)
         d(dzengjiahe3(i,1),dzengjiahe3(i,2))=inf;
         d(dzengjiahe3(i,2),dzengjiahe3(i,1))=inf;
     end
  end
  
   for h=1:length(jilu)
       if abs(d(jilu(h,1),jilu(h,2)))<zhibiao 
           d(jilu(h,1),jilu(h,2))=1;
           d(jilu(h,2),jilu(h,1))=1;
       else
           d(jilu(h,1),jilu(h,2))=inf;  
           d(jilu(h,2),jilu(h,1))=inf;
           dzengjia=[dzengjia;jilu(h,1),jilu(h,2);jilu(h,2),jilu(h,1)];   
       end
   end
   if j==1
      dzengjiahe1=[dzengjiahe1;dzengjia];
   elseif j==2
      dzengjiahe2=[dzengjiahe2;dzengjia];
   elseif j==3
       dzengjiahe3=[dzengjiahe3;dzengjia];
   end
   d(find(d~=1&d~=inf))=1;
   dhalf=sqrt((lon-lon').^2+(lat-lat').^2);
   d3=d.*dhalf;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tree_exist=0;
dc=0.1;
while tree_exist==0
    dhalf(dhalf>dc)=0;
    dhalf=sparse(tril(dhalf));
    [Tree] = graphminspantree(dhalf);
    [index1_path,index2_path]=find(Tree);
    if isempty(index1_path)==1||length(index1_path)<n-1
        dc=dc*2;
    else
        tree_exist=1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tree_exist=0;
dc=0.1;
while tree_exist==0
    d1=sparse(tril(d1));
    [Tree] = graphminspantree(d1);
    [index1_path1,index2_path1]=find(Tree);
    if isempty(index1_path1)==1||length(index1_path1)<n-1
        dc=dc*2;
    else
        tree_exist=1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tree_exist=0;
dc=0.1;
while tree_exist==0    
    d2=sparse(tril(d2));
    [Tree] = graphminspantree(d2);
    [index1_path2,index2_path2]=find(Tree);
    if isempty(index1_path2)==1||length(index1_path2)<n-1
        dc=dc*2;
    else
        tree_exist=1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tree_exist=0;
dc=0.1;
while tree_exist==0
    d3=sparse(tril(d3));
    [Tree] = graphminspantree(d3);
    [index1_path3,index2_path3]=find(Tree);
    if isempty(index1_path3)==1||length(index1_path3)<n-1
        dc=dc*2;
    else
        tree_exist=1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1
bianjie=[];
for i=1:length(index1_path)
       if abs(beta(index1_path(i),1)-beta(index2_path(i),1))>1e-3
           bianjie=[bianjie;index1_path(i),index2_path(i)];
       end
end
count_y1=length(bianjie);%%%%%%%%计算边界处的边数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2
bianjie=[];
for i=1:length(index1_path)
       if abs(beta(index1_path(i),2)-beta(index2_path(i),2))>1e-3
           bianjie=[bianjie;index1_path(i),index2_path(i)];
       end
end
count_y2=length(bianjie);%%%%%%%%计算边界处的边数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3
bianjie=[];
for i=1:length(index1_path)
       if abs(beta(index1_path(i),3)-beta(index2_path(i),3))>1e-3
           bianjie=[bianjie;index1_path(i),index2_path(i)];
       end
end
count_y3=length(bianjie);%%%%%%%%计算边界处的边数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4
bianjie=[];
for i=1:length(index1_path1)
       if abs(beta(index1_path1(i),1)-beta(index2_path1(i),1))>1e-3
           bianjie=[bianjie;index1_path1(i),index2_path1(i)];
       end
end
count_1=length(bianjie);%%%%%%%%计算边界处的边数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
bianjie=[];
for i=1:length(index1_path2)
       if abs(beta(index1_path2(i),2)-beta(index2_path2(i),2))>1e-3
           bianjie=[bianjie;index1_path2(i),index2_path2(i)];
       end
end
count_2=length(bianjie);%%%%%%%%计算边界处的边数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%6    
 bianjie=[];
for i=1:length(index1_path3)
       if abs(beta(index1_path3(i),3)-beta(index2_path3(i),3))>1e-3
           bianjie=[bianjie;index1_path3(i),index2_path3(i)];
       end
end
count_3=length(bianjie);%%%%%%%%计算边界处的边数
count_yuanshi=[count_y1,count_y2,count_y3];
count_hou=[count_1,count_2,count_3];
count_yuanshi_z=sum(count_yuanshi);
count_hou_z=sum(count_hou);
end