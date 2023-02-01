function[H,Hhou,dzengjiahe1hou,dzengjiahe2hou,dzengjiahe3hou]=SCC_spanning_treegai2fenpro(lon,lat,p,dc,x,y,lamuda,beta_hat_SCC1,t,Hqian,dzengjiahe1,dzengjiahe2,dzengjiahe3,zhibiaocishu)
% Compute n*(n-1) matrix constructed from edge set of minimum spanning tree
% INPUT ARGUMNETS:
% lon: [n,1] matrix storing the longitudes of n locations
% lat: [n,1] matrix storing the latitudes of n locations
% p  : number of variables
% dc : a critical distance larger than which the two vertices will not be
%      connected. This value is used only for the reduction of computational
%      burden.

%OUTPUT ARGUMENTS
% H  :  n*(n-1) matrix constructed from edge set of minimum spanning tree

n=length(lon);

dbeta=squeeze(beta_hat_SCC1(t,:,:));%%%%%%%%%%%%%%%%%%%%%%%%

dhalf=sqrt(((lon-lon').^2+(lat-lat').^2));

H=zeros((n-1)*p,n*p);

for j=1:p
    j;


tree_exist=0;
    %v1=normalize(x(:,1));
    %v2=normalize(x(:,2));
    %v3=normalize(x(:,3));%%%全1
    %v4=normalize(y(:)); 
   
   %[m,n]=size(dhalf);
   %for g=1:m
       %a=sum(exp(dhalf(g,:)));
       %for h=1:n
           %dhalf(g,h)=exp(dhalf(g,h))/a;
       %end
   %end
   
   
   %dhalf=abs(lon-lon')+abs(lat-lat');
   %maxd=max(max(dhalf));
   %dhalf=dhalf./maxd;
   %dhalf=exp(dhalf);

   dbeta(:,j)=dbeta(:,j)/std(dbeta(:,j));%%%%%%%%%%%%%%%%%%%%%
   d=abs((dbeta(:,j)-dbeta(:,j)'));%%%%%%%%%%%%%%%%%%%%%%%%
   range=max(max(d))-min(min(d));
   shangjie=999999;

   %a=find(d>zhibiao);
   %length(a)
   [m,n]=size(d);   
   
   

   Hqianbianliang=Hqian((n-1)*(j-1)+1:(n-1)*j,n*(j-1)+1:n*j);
   jilu=[];
   for l=1:n-1
       jilu=[jilu;find(Hqianbianliang(l,:)~=0)];
   end
jilu;
  % for i=1:m
     %  for t=1:n
        %  if d(i,t)<=zhibiao(i,i)
           %   d(i,t)=1;
         % else
             % d(i,t)=shangjie;
         % end
      % end
  % end
  dzengjia=[];%%%%%%%%%%%%
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
       if abs(d(jilu(h,1),jilu(h,2)))<zhibiao  %%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%
           d(jilu(h,1),jilu(h,2))=1;%%%%%%%%%%%%%%%%%%%%
           d(jilu(h,2),jilu(h,1))=1;%%%%%%%%%%%%%%%%%%%%
       else
           d(jilu(h,1),jilu(h,2))=inf;  %%%%%%%%%%%%%%%%%%%%%永久断掉的
           d(jilu(h,2),jilu(h,1))=inf;
           dzengjia=[dzengjia;jilu(h,1),jilu(h,2);jilu(h,2),jilu(h,1)];    %%%%%%%%%%%%
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
   %d(find(d~=shangjie&d~=1&d~=inf))=shangjie;
   %huanyuan=find(d==shangjie);
   %d(find(d==shangjie))=shangjie;
   d=d.*dhalf;%%%%%
   %d(huanyuan)=shangjie;
   %
   %dcanshu=abs((dbeta(:,j)-dbeta(:,j)'));
   %jingxi=find(d==1);
   %d(jingxi)=dcanshu(jingxi);
  d;
length(find(d==inf))%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   
while tree_exist==0
    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%


      
       %for i=1:m
           %d(i,:)=d(i,:)/min(d(i,:));
       %end
       
       
    %range=max(max(d))-min(min(d));
    %d(find(d==0))=1;
    %d(find(d~=0))=9999;
    %d=exp(d/range);
    %d=abs(dbeta(:,j)-dbeta(:,j)');
    %d=abs(dbeta(:,1)-dbeta(:,1)')+abs(dbeta(:,2)-dbeta(:,2)')+abs(dbeta(:,3)-dbeta(:,3)');
    %d=sqrt((dbeta(:,1)-dbeta(:,1)').^2+(dbeta(:,2)-dbeta(:,2)').^2+(dbeta(:,3)-dbeta(:,3)').^2);
    %d=dhalf;
 
    %(1-lamuda)*((v1-v1').^2+(v2-v2').^2+(v4-v4').^2);
    %(v1-v1').^2+(v2-v2').^2+(v4-v4').^2;
    
    %[q,r]=size(d);
    %d2=sqrt((lon-lon').^2+(lat-lat').^2);
    %pandingzhi=1*mean(mean(d2));
    %d(d2>pandingzhi)=inf;
    
    
   % d(d>dc)=0;
    
    d=sparse(tril(d));
    
    [Tree] = graphminspantree(d);
    [index1_path,index2_path]=find(Tree);
    if isempty(index1_path)==1||length(index1_path)<n-1
        dc=dc*2;
    else
        tree_exist=1;
    end
end

   % plot(lon,lat,'o')
   % hold on 
   % for i=1:length(index1_path)
        %    x111=lon(index1_path(i,1),1);
         %   y111=lat(index1_path(i,1),1);
         %   x222=lon(index2_path(i,1),1);
         %   y222=lat(index2_path(i,1),1);
         %   x=[x111,x222];%待连接的俩个点的横坐标
         %   y=[y111,y222];%待连接的俩个点的坐标
         %   mst=line(x,y);
         %   mst.Color='red';
          %  hold on
  %  end




H1=zeros(n-1,n);

ppp=sub2ind(size(H1),[1:n-1]',index1_path);
H1(ppp)=1;
ppp=sub2ind(size(H1),[1:n-1]',index2_path);
H1(ppp)=-1;



%%%%%%%%循环位置
    H((j-1)*(n-1)+1:j*(n-1),(j-1)*n+1:j*n)=H1;
end 

for j=1:p
    H(p*(n-1)+j,1+(j-1)*n:j*n)=1/n;
end

dzengjiahe1hou=dzengjiahe1;

dzengjiahe2hou=dzengjiahe2;
dzengjiahe3hou=dzengjiahe3;

Hhou=H;
end