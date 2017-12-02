function [ Init_leaks ] = Initializ( cnt,n,m )
% %Author: Mozhou Gao
% %University of Calgary 
% %MGIS student in Department of Geography
% %Input: number of wells in each cell, dimension of grids
% %Output: The 3D data box 
%% generate the initial leaks 
mu = 6.*cnt; %caculate the mean for each cell
Init_leaks = zeros(n,m);
for i =1:n
    for j = 1:m 
        tempL = 0:1:(12*cnt(i,j));%form the distribution based on cell density
        sig = std(tempL);%calculate the standard deviation  
        Init_leaks (i,j) = round(normrnd(mu(i,j),sig));%generate the leaks based on the distribution of each cell
        if Init_leaks(i,j) < 0 
           Init_leaks(i,j) = 0; %truncate the lower boundary
        end
        
    end 
end 
% error check
checka = (cnt >0);
checkb = (Init_leaks > 0);
checkc = (cnt == 0 &  Init_leaks ~= 0);
erLocation = NaN(2,1);
for i = 1:n
    for j = 1:m 
        if checkc(i,j) ~= 0
            sprintf('Error: The Random Leaks from nowhere');
            erLocation (1) = i;
            erLocation (2) = j;
        end 
    end 
end 
end

