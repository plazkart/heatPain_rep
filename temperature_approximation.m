function y_1 = temperature_approximation(temp, VAS, selected_temp)
% 
% temp = [39, 42, 45, 46, 47, 48];
% VAS = [2, 3, 4, 5, 6, 7];


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

%fitting with polynom fucntion
% N=2;
% p = polyfit(VAS, temp, N);
% x1 = linspace(min(VAS), max(VAS), 100);
% y1 = polyval(p,x1);
% fitting with sigmoid func
temp(end+1) = 32; VAS(end+1) = 0;
K = 10; A= 0;
fitfun = fittype( @(X50, B, x) (K-A)./(1+exp(-B*(x-X50))));
% fitted_model = fit(temp',VAS',fitfun, 'StartPoint', [1, 42]);
fitted_model = fit(temp',VAS',fitfun, 'Lower', [40, 0], 'Upper', [51, 1.5]);
x1 = linspace(min(temp), max(temp), 100);
y1 = fitted_model(x1);
% coeffvals = coeffvalues(fitted_model);

figure(); plot(VAS, temp, 'or', x1, y1, '-b');
grid on
xlim([30, 50]);
ylim([0, 10]);

y_1 = fitted_model(selected_temp);

% x_val = 7;
% % y_pred = p(1)*(x_val)^2 + p(2)*x_val + p(3)
% syms x
% eqn = p(1)*x^2 + p(2)*x + p(3) == selected_temp;
% eqn = p(1) == selected_temp;
% y_1 = solve(eqn, x);
% y_1 = double(y_1);
% y_1 = y_1(1);

