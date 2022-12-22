% error contour plot
load('Error_matrix_train.mat', 'error_matrix_train')
load('Error_matrix_test.mat', 'error_matrix_test')
alpha_Gin_vals = [1, 2, 4, 8, 16];
alpha_Gres_vals = [0.125, 0.25, 0.5, 1, 2, 4];
[X, Y] = meshgrid(alpha_Gres_vals, alpha_Gin_vals)
contour(X, Y, error_matrix_test) 
xlabel('$\alpha_{Gres}$','Interpreter',"latex", 'FontSize',30)
ylabel('$\alpha_{Gin}$','Interpreter',"latex", 'FontSize',30)

%% input raster plots
data_group = cell(10);
for i = 0:9
    data_group{i+1} = DATA([DATA.type] == i);
end


digit = 9; % replace with slider in livescript (0 to 9)
size_digit = size(data_group{digit+1}, 2);
data_digit = data_group{digit+1};

[a,b] = find(data_digit(1).S);
figure(2);plot(b,a,'.')

[a,b] = find(data_digit(end-2).S);
figure(2);plot(b,a,'.')

%% output raster plots
data_group = cell(10);
for i = 0:9
    data_group{i+1} = DATA([DATA.type] == i);
end

digit = 9; % replace with slider in livescript (0 to 9)
size_digit = size(data_group{digit+1}, 2);
data_digit = data_group{digit+1};

[a,b] = find(data_digit(k).RES);
figure(2);
% plot(mean(data_digit(k).RES, 2))
plot(b,a,'.')