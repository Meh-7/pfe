data = table2array(readtable('databp2.csv', 'HeaderLines', 0));
theta = mean_theta_hat;
f = likelihooddelt(theta,data);
x = -2:1e-08:2;
y = f(x);

