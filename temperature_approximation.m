temp = [39, 42, 45, 46, 47, 48];
VAS = [2, 3, 4, 5, 6, 7];


% plot(VAS, temp, 'or');
% grid on
% xlim([0, 10]);
% ylim([30, 50]);

% N=2;
% p = polyfit(temp, VAS, N);
% x1 = linspace(min(temp), max(temp), 100);
% y1 = polyval(p,x1);
% plot(temp, VAS, 'or', x1, y1, '-');
% grid on
% xlim([30, 50]);
% ylim([0, 10]);

% x(y==6)
% x1(y1==6)

N=2;
p = polyfit(VAS, temp, N);
x1 = linspace(min(VAS), max(VAS), 100);
y1 = polyval(p,x1);
plot(VAS, temp, 'or', x1, y1, '-b');
grid on
xlim([0, 10]);
ylim([30, 50]);

x_val = 7;
y_pred = p(1)*(x_val)^2 + p(2)*x_val + p(3)


