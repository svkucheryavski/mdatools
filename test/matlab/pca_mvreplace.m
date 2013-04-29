clear

%% simple example
x = 1:10
data = [x' 2*x' 3*x'];
data(5, 1) = NaN;
data(8, 3) = NaN;

options = mdcheck('options');
options.display = 'on';

[a, b, newdata] = mdcheck(data, options)

newdata

%% People example

%load('People.mat')