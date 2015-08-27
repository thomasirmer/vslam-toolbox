% load data
load('Obs.mat');

% homogenous line coordinates
l = Obs.exp.e(1);
m = Obs.exp.e(2);
n = Obs.exp.e(3);

% line segment coordinates
line_seg_x = [Obs.meas.y(1), Obs.meas.y(3)];
line_seg_y = [Obs.meas.y(2), Obs.meas.y(4)];

% homogenous line
line_homo_x = 0:600;
line_homo_y = - (n + l * line_homo_x) / m;

% standard deviation for each coordinate
s_dev_l = 0.25; % 3 * sqrt(Obs.exp.E(1,1));
s_dev_m = 0.05; % 3 * sqrt(Obs.exp.E(2,2));
s_dev_n = 10; % 3 * sqrt(Obs.exp.E(3,3));

% lines with standard deviation
line_homo_dev_l_pos = - (n + (l + s_dev_l) * line_homo_x) / m;
line_homo_dev_l_neg = - (n + (l - s_dev_l) * line_homo_x) / m;

line_homo_dev_m_pos = - (n + l * line_homo_x) / (m + s_dev_m);
line_homo_dev_m_neg = - (n + l * line_homo_x) / (m - s_dev_m);

line_homo_dev_n_pos = - ((n + s_dev_n) + l * line_homo_x) / m;
line_homo_dev_n_neg = - ((n - s_dev_n) + l * line_homo_x) / m;

% plot
hold off
plot(line_seg_x, line_seg_y, 'black*--');
hold on;
plot(line_homo_x, line_homo_y, 'black--');

plot(line_homo_x, line_homo_dev_l_pos, 'red');
plot(line_homo_x, line_homo_dev_l_neg, 'red');

plot(line_homo_x, line_homo_dev_m_pos, 'green');
plot(line_homo_x, line_homo_dev_m_neg, 'green');

plot(line_homo_x, line_homo_dev_n_pos, 'blue');
plot(line_homo_x, line_homo_dev_n_neg, 'blue');

axis([150, 300, 100, 300]);
legend('line segment', ...
       'homogenous line', ...
       'std. dev. l pos',  ...
       'std. dev. l neg', ...
       'std. dev. m pos', ...
       'std. dev. m neg', ...
       'std. dev. n pos', ...
       'std. dev. n neg');
