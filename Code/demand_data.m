function out = demand_data(P_dem_max)

index = linspace(1,23,23);
rem = 6; %buildings to remove based on negative energy consumption
A = setdiff(index,rem);
new_con = [];
for i=1:length(A)
    rng = xlsColNum2Str(A(i)+3); %start at column D
    searchCol = strcat(rng,'2:',rng,'8569'); %search data column from row 2 to 8569
    dataCol = xlsread('HW4_Train_Data.csv','HW4_Train_Data',searchCol{1}); ...
        %imported data column
    new_con= [new_con, dataCol./max(dataCol)]; %normalize data to max energy consumption
end
num_buildings = length(index)-length(rem);
B = permute(reshape(new_con, [24 7 51 num_buildings]), [4 3 2 1]); ...
    %num_buildings, 51 weeks of data, 7 days/week, 24 hours/day
hourly_avg = zeros(24,1);
for i=1:7
    for j=1:num_buildings
        for k=1:51
            for l=1:24
                hourly_avg(l) = mean(mean(B(:,:,i,l)));
            end
        end
    end
end

% Compute demand output
out = P_dem_max*hourly_avg;