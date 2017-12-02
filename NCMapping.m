
%call WeatherMapping Function 
[ freq,Tr,Wr,lonData,latData,timeData ] = WeatherMapping( 3,7,0.5 );
lo = length(lonData);
la = length(latData);
day = [ 365 365 365 366 365 365 365 366 365 365 365 366 365 365 365 366];

OGIday = zeros(lo,la,length(day));
MWGsday = zeros(lo,la,length(day));
index1 = 1;
index2 = 365;
for k = 1:1:16
    YearOGI = freq(:,:,index1:index2);
    OGIday(:,:,k) = sum(YearOGI,3);
    for i = 1:1:lo
        for j = 1:1:la
            %find the index of evary OGI working day 
             OGIdayIndex = find(YearOGI(i,j,:));
             L =length(OGIdayIndex);%check the number of OGI days for each cell
             if L == 0 %if OGI days = 0 
                MWGsday(i,j,k)=365 ;
             elseif L ==1 %if OGI days = 1
                 tempM = YearOGI(i,j,:); 
                MWGsday(i,j,k)=365 - OGIdayIndex; %Max Winter gap = 366
             else %otherwise
                 Diff = diff(OGIdayIndex);% calculate the differences between every OGI days index 
                 FD = OGIdayIndex(1);
                 LD = 365 - OGIdayIndex(L);
                 MaxGap = max(Diff);%Find the maximum difference
                 TempList = [FD LD MaxGap];
                 MWGsday(i,j,k) = max(TempList);
             end 
        end
    end
    index1 = index1 + day(k);
    index2 = index2 + day(k);
end
AveOGI = (sum(OGIday,3)./16);
StdOGI = std(OGIday,0,3);
AveMWG = (sum(MWGsday,3)./16);
StdMWG = std(MWGsday,0,3);

%% Export OGIDays
%Average
nccreate('OGID_mean_T3_W7.nc','OGIdays','Dimensions',{'lon',lo,'lat',la},'Format','classic');
nccreate('OGID_mean_T3_W7.nc','lon','Dimensions',{'lon',lo},'Format','classic');
nccreate('OGID_mean_T3_W7.nc','lat','Dimensions',{'lat',la},'Format','classic');
ncwrite('OGID_mean_T3_W7.nc','OGIdays',AveOGI);
ncwrite('OGID_mean_T3_W7.nc','lon',lonData);
ncwrite('OGID_mean_T3_W7.nc','lat',latData);
%Standard Deviation
nccreate('OGID_std_T3_W7.nc','stdOGI','Dimensions',{'lon',lo,'lat',la},'Format','classic');
nccreate('OGID_std_T3_W7.nc','lon','Dimensions',{'lon',lo},'Format','classic');
nccreate('OGID_std_T3_W7.nc','lat','Dimensions',{'lat',la},'Format','classic');
ncwrite('OGID_std_T3_W7.nc','stdOGI',StdOGI);
ncwrite('OGID_std_T3_W7.nc','lon',lonData);
ncwrite('OGID_std_T3_W7.nc','lat',latData);

%% Export MWGs
nccreate('MWG_mean_T3_W7.nc','GapLength','Dimensions',{'lon',lo,'lat',la},'Format','classic');
nccreate('MWG_mean_T3_W7.nc','lon','Dimensions',{'lon',lo},'Format','classic');
nccreate('MWG_mean_T3_W7.nc','lat','Dimensions',{'lat',la},'Format','classic');
ncwrite('MWG_mean_T3_W7.nc','GapLength',AveMWG);
ncwrite('MWG_mean_T3_W7.nc','lon',lonData);
ncwrite('MWG_mean_T3_W7.nc','lat',latData);
% 
%standard Deviation
nccreate('MWG_std_T3_W7.nc','stdMWG','Dimensions',{'lon',lo,'lat',la},'Format','classic');
nccreate('MWG_std_T3_W7.nc','lon','Dimensions',{'lon',lo},'Format','classic');
nccreate('MWG_std_T3_W7.nc','lat','Dimensions',{'lat',la},'Format','classic');
ncwrite('MWG_std_T3_W7.nc','stdMWG',StdMWG);
ncwrite('MWG_std_T3_W7.nc','lon',lonData);
ncwrite('MWG_std_T3_W7.nc','lat',latData);
fprintf('Data Exporting Finish!!')
