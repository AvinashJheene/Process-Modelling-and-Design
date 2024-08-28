clear

%% Ethanol - water vapor-liquid equilibrium %%

x = [0 0.0186 0.0476 0.0673 0.0881 0.1102 0.1424 0.1894 0.2069]; % mole fraction 
y = [0 0.0105 0.0272 0.0375 0.0492 0.0624 0.0809 0.1078 0.1182]; % mole fraction
X = x./(1 - x);                                                 % mole ratio
Y = y./(1 - y);                                                 % mole ratio

subplot(3,2,1);
scatter(X,Y);
xlabel('X');
ylabel('Y');
hold on;

%% VLE fitting %%

F_fit = @(x_fit,x_data) (x_fit(1).*x_data)./(1 + (x_fit(2).*x_data));       % Y-X fitting
x_fit0 = [1 1];
[x_fit, resnorm] = lsqcurvefit(F_fit,x_fit0,X,Y);   % Fitting method 1

beta = nlinfit(X,Y,F_fit,x_fit0);                   % Fitting method 2

% Fitted constants a = 0.5576, b = 0.3287

X_eq = linspace(0,0.50,1000);

% Fitted curve Y_eq
Y_eq = (x_fit(1).*X_eq)./(1 + (x_fit(2).*X_eq));   
subplot(3,2,1);
plot(X_eq,Y_eq,Color='b');                                               % Use Y-eq from now on
xlabel('X_{eq}');
ylabel('Y_{eq}');
hold on; 

%% Gas feed composition and solute reomval %%
% Gas feed rate - 2000 kg/hr
% CO2 - 85% mol, EtOH - 15% mol

Mav = 0.85*44 + 0.15*46;                                        % Avg molecular weight of feed gas = 44.3
G1 = 2000/Mav;                                                  % gas feed rate kmol/hr = 45.1467
%G1 = 43.76;                                                    % If someone has obtained G1 = 43.76 that is also okay
y1 = 0.15;                                                      % feed concentration (mol fraction), y1 = 0.15
Y1 = y1/(1-y1);                                                 % feed concentration (mol ratio), Y1 = 0.1765
Gs = G1*(1-y1);                                                 % feed rate solute-free basis, Gs = 38.3747 kmol/hr
G1_etoh = G1*y1;                                                % etoh entering, G1_etoh = 6.7720 kmol/hr

solute_removal = [0.86 0.88 0.92 0.94 0.96 0.98];               % array to define ethanol recoveries
% The first two values of the array are for the bonus question

%% Solvent composition %%
% EtOH in solvent feed - 0.0 mole%

x2 = 0.0;                                                       % x2 and X2 will be the same for all
X2 = x2/(1-x2);

% Full marks if x2 and X2 has been defined in the code atleast once
%% Array initialization %%

sz = size(solute_removal,2);
G2_etoh = zeros(1,sz);
Y2 = zeros(1,sz);
y2 = zeros(1,sz);
X1_max = zeros(1,sz);
Ls_min = zeros(1,sz);
Ls = zeros(1,sz);
NTU_1 = zeros(1,sz);
NTU_2 = zeros(1,sz);
pinchSlope = zeros(1,sz);
pinchpoint = zeros(2,sz);
Ht = zeros(1,sz);

%% Exit gas concentration %%

for i = 1:sz
    G2_etoh(i) = G1_etoh*(1 - solute_removal(i));                                         % etoh leaving
    Y2(i) = G2_etoh(i)/Gs;                                                                % exit etoh conc.
    y2(i) = Y2(i)/(1+Y2(i));                                                              % exit etoh conc.
end

% G2_etoh (kmol/hr) : 92% - 0.5418, 94% - 0.4063, 96% - 0.2709, 98% - 0.1354
% Y2: 92% - 0.0141, 94% - 0.0106, 96% - 0.0071, 98% - 0.0035
% y2: 92% - 0.0139, 94% - 0.0105, 96% - 0.0070, 98% - 0.0035
% x2: 0 for all
% X2: 0 for all
% For bonus question: Y2: 86% - 0.0247, 88% - 0.0212
% y2: 86% - 0.0241, 88% - 0.0207

%% Minimum and actual solvent flow rate %%

pinch0 = [0.1 0.1];                                                                     % initial guess for determining pinch point
for i = 1:sz              
    f_minSolvent = @(ct) minSolvent(ct, x_fit(1), x_fit(2), X2, Y2(i));                 % calling the function to determine pinch point
    pinch = fsolve(f_minSolvent,pinch0);
    if (pinch(2) >= Y1)
        pinch(2) = Y1;
        pinch(1) = Y1/(x_fit(1) - x_fit(2)*Y1);
    end
    pinchpoint(1,i) = pinch(1);
    pinchpoint(2,i) = pinch(2);
    pinchSlope(i) = (pinch(2) - Y2(i))/(pinch(1) - X2);
    X1_max(i) = ((Y1 - Y2(i))/pinchSlope(i)) + X2;                                      
    Ls_min(i) = Gs*pinchSlope(i);                                                         % Minimum solvent rate
    Ls(i) = 1.25*Ls_min(i);                                                               % Actual solvent rate
end
% Pinch point coordinates:
% 92%: x-0.3054, y-0.1548; 
% 94%: x-0.2610, y-0.1340; 
% 96%: x-0.2098, y-0.1094; 
% 98%: x-0.1454, y-0.0774;

% Minimumn solvent rate, Ls_min:
% 92%: 17.6720;			 
% 94%: 18.1505; 
% 96%: 18.7264; 
% 98%: 19.4903;

subplot(3,2,2);
plot(solute_removal,Ls_min);
xlabel('recovery');
ylabel('L_{s,min}');

% Minimum solvent rate increases with the required recovery as more solvent
% needs to be provided for higher solute removals. All points including
% bonus question recoveries are included in this plot.

% For bonus question:
% Pinch point coordinates:
% 86%: x-0.3532, y-0.3532; 
% 88%: x-0.1765, y-0.1765;

% Minimumn solvent rate, Ls_min:
% 86%: 16.4880;			 
% 88%: 16.8715;
%% Operating line %%

operatingLine = zeros(size(X_eq,2),sz);
minSolvent_operatingLine = zeros(size(X_eq,2),sz);

for i=1:sz
    for j=1:size(X_eq,2)
        operatingLine(i,j) = Y2(i) + (Ls(i)/Gs)*(X_eq(j) - X2);                         % Operating line for all recoveries
        minSolvent_operatingLine(i,j) = Y2(i) + (Ls_min(i)/Gs)*(X_eq(j) - X2);          % Minimum solvent operating lines for all recoveries
    end
end

% Give full marks for part (vii) if atleast one operating line is shown in the plot

subplot(3,2,3)
plot(X_eq,Y_eq,Color='b'); 
hold on
plot(X_eq, operatingLine(5,:), "Color","r","LineStyle","--");               % Operating line for 94% recovery
hold on
yline(Y1,Color='g',LineStyle='-.');
hold on
plot(X_eq, minSolvent_operatingLine(5,:));                                  % Minimum solvent operating line for 94% recovery
axis([0 0.7 0 0.4]);
hold on
plot(pinchpoint(1,4),pinchpoint(2,4), Marker="*")
xlabel('X');
ylabel('Y');
legend('Equilibrium Curve','Operating line','Y1', 'Min solv op line', 'Pinch Point')

X1_max_94 = X1_max(4);

%% Bonus question

% For 86% recovery

subplot(3,2,4)
plot(X_eq,Y_eq,Color='b'); 
hold on
plot(X_eq, operatingLine(1,:), "Color","r","LineStyle","--");               % Operating line for 86% recovery
hold on
yline(Y1,Color='g',LineStyle='-.');
hold on
plot(X_eq, minSolvent_operatingLine(1,:));                                  % Minimum solvent operating line for 86% recovery
axis([0 0.7 0 0.4]);
hold on
plot(pinchpoint(1,1),pinchpoint(2,1), Marker="*")
xlabel('X');
ylabel('Y');
legend('Equilibrium Curve','Operating line 86% rec','Y1','Min solv op line', 'Pinch Point')


% For 88% recovery

subplot(3,2,5)
plot(X_eq,Y_eq,Color='b'); 
hold on
plot(X_eq, operatingLine(2,:), "Color","r","LineStyle","--");               % Operating line for 88% recovery
hold on
yline(Y1,Color='g',LineStyle='-.');
hold on
plot(X_eq, minSolvent_operatingLine(2,:));                                  % Minimum solvent operating line for 86% recovery
axis([0 0.7 0 0.4]);
hold on
plot(pinchpoint(1,2),pinchpoint(2,2), Marker="*")
xlabel('X');
ylabel('Y');
legend('Equilibrium Curve','Operating line 88% rec','Y1','Min solv op line', 'Pinch Point')