all_x_data = load('all_x.mat'); 
all_x = all_x_data.all_x;
all_y_data = load('all_y.mat'); 
all_y = all_y_data.all_y;
SA_x_data = load('SA_x.mat'); 
SA_x = SA_x_data.SA_x;
SA_y_data = load('SA_y.mat'); 
SA_y = SA_y_data.SA_y;

for i=1:100
    SA_x(i) = 1000*(i-1)+1;
end

plot(all_x,all_y,'k.',SA_x,SA_y,'c*')
