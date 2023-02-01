function[H,Hhou,dzengjiahe1hou,dzengjiahe2hou,dzengjiahe3hou]=SCC_spanning_treegai(lon,lat,p,dc,x,y)
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


tree_exist=0;
  
    v1=normalize(x(:,1));
    v2=normalize(x(:,2));
    v3=normalize(x(:,3));%%%全1
    v4=normalize(y(:));
while tree_exist==0
    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%

  
    d=sqrt((lon-lon').^2+(lat-lat').^2);
    %d=sqrt(lamuda*((lon-lon').^2+(lat-lat').^2)+(1-lamuda)*((v1-v1').^2+(v2-v2').^2+(v4-v4').^2));
    d(d>dc)=0;
    %(1-lamuda)*((v1-v1').^2+(v2-v2').^2+(v4-v4').^2);
    %(v1-v1').^2+(v2-v2').^2+(v4-v4').^2;
    
    
    d=sparse(tril(d));
    
    [Tree] = graphminspantree(d);
    [index1_path,index2_path]=find(Tree);
    if isempty(index1_path)==1||length(index1_path)<n-1
        dc=dc*2;
    else
        tree_exist=1;
    end
end
   
  %  plot(lon,lat,'.','MarkerSize',10,'LineWidth',5)
  %  hold on 
   % for i=1:length(index1_path)
    %        x111=lon(index1_path(i,1),1);
     %       y111=lat(index1_path(i,1),1);
      %      x222=lon(index2_path(i,1),1);
      %      y222=lat(index2_path(i,1),1);
       %     x=[x111,x222];%待连接的俩个点的横坐标
       %     y=[y111,y222];%待连接的俩个点的坐标
       %     mst=line(x,y);
        %    mst.Color='red';
       %     hold on
   % end




H1=zeros(n-1,n);

ppp=sub2ind(size(H1),[1:n-1]',index1_path);
H1(ppp)=1;
ppp=sub2ind(size(H1),[1:n-1]',index2_path);
H1(ppp)=-1;

H=zeros((n-1)*p,n*p);
for j=1:p
    H((j-1)*(n-1)+1:j*(n-1),(j-1)*n+1:j*n)=H1;
end 

for j=1:p
    H(p*(n-1)+j,1+(j-1)*n:j*n)=1/n;
end
    