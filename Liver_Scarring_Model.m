% Written by James Guentert

clc; clear; close all;

%% Parameters (using L and L/min consistently)

p.Qco   = 5;            % L/min (cardiac output)
p.QH    = 0.25*p.Qco;   % L/min hepatic blood flow
p.VL    = 1.7;          % L liver volume
p.VB    = 20 - p.VL;    % L effective blood/water volume
p.VD    = 1.05;         % Liver density kg/L

% Kinetics (M/min and M)
p.VmaxE = 3.256E-3 *1/p.VL;          % M/min
p.KmE   = 0.8183E-3;                 % M
p.VmaxA = mean([2.4,4.7])*p.VD*1E-3; % M/min
p.KmA   = 1E-3;                      % M

%% Input u(t): mol/min into blood

tOn = 60;                   % drink for one hour (in minutes)
tEnd = tOn*6;              % analyze for one day (in minutes)
uin_mol_per_min = 0.1*0.0152;   % mol/min, 3 standard drink/hour
tU  = 0:1:tEnd;             % time vector
ut  = uin_mol_per_min*(tU <= tOn); % input u(t)

% Variety of alpha values
alphas = [1, 0.9, 0.75, 0.5, 0.3];

% Cooresponding Colors
colors={'#1a1880','#18801a','#FCE303','#CF5B08','#BD0808'};


%% Ploting

figure('Color',hex2rgb('#212121'));

% E_B overlay
subplot(2,2,1); hold on; grid on
set(gca,'Color',hex2rgb('#adadad'))
title('Blood Ethanol E_B(t)','Color','w')
xlabel('Time (hr)'); ylabel('Concentration [mM]')

% E_L overlay
subplot(2,2,2); hold on; grid on
set(gca,'Color',hex2rgb('#adadad'))
title('Liver Ethanol E_L(t)','Color','w')
xlabel('Time (hr)'); ylabel('Concentration [mM]')

% A_B overlay
subplot(2,2,3); hold on; grid on
set(gca,'Color',hex2rgb('#adadad'))
title('Blood Acetaldehyde A_B(t)','Color','w')
xlabel('Time (hr)'); ylabel('Concentration [mM]')

% A_L overlay
subplot(2,2,4); hold on; grid on
set(gca,'Color',hex2rgb('#adadad'))
title('Liver Acetaldehyde A_L(t)','Color','w')
xlabel('Time (hr)'); ylabel('Concentration [mM]')

for k = 1:numel(alphas)
    p.alpha = alphas(k);
    p.color = colors(k);

    % y = [E_B; E_L; A_B; A_L] in M
    y0 = [0; 0; 0; 0];

    odefun = @(t,y) rhs(t,y,p,tU,ut);
    [t,y] = ode45(odefun, [0 tEnd], y0);

    E_B = y(:,1);
    E_L = y(:,2);
    A_B = y(:,3);
    A_L = y(:,4);
    
    % Save
    sol(k).alpha = p.alpha;
    sol(k).t = t;
    sol(k).y = y;

    % Plot (and convert M -> mM)
    subplot(2,2,1);
    plot(t/60, E_B*1e3, 'LineWidth', 2, 'DisplayName', sprintf('\\alpha = %.1f', p.alpha),'Color',hex2rgb(p.color));

    subplot(2,2,2);
    plot(t/60, E_L*1e3, 'LineWidth', 2, 'DisplayName', sprintf('\\alpha = %.1f', p.alpha),'Color',hex2rgb(p.color));
    
    subplot(2,2,3);
    plot(t/60, A_B*1e3, 'LineWidth', 2, 'DisplayName', sprintf('\\alpha = %.1f', p.alpha),'Color',hex2rgb(p.color));

    subplot(2,2,4);
    plot(t/60, A_L*1e3, 'LineWidth', 2, 'DisplayName', sprintf('\\alpha = %.1f', p.alpha),'Color',hex2rgb(p.color));

end

subplot(2,2,1); legend('Location','best'); xticks(0:2:tEnd); yticks(0:1:180);
ax = gca; ax.XColor = 'w'; ax.YColor = 'w'; ylim([0 inf])
subplot(2,2,2); legend('Location','best'); xticks(0:2:tEnd); yticks(0:1:180);
ax = gca; ax.XColor = 'w'; ax.YColor = 'w'; ylim([0 inf])
subplot(2,2,3); legend('Location','best'); xticks(0:2:tEnd); yticks(0:0.125:180);
ax = gca; ax.XColor = 'w'; ax.YColor = 'w'; ylim([0 inf])
subplot(2,2,4); legend('Location','best'); xticks(0:2:tEnd); yticks(0:0.25:180);
ax = gca; ax.XColor = 'w'; ax.YColor = 'w'; ylim([0 inf])

% RHS function
function dydt = rhs(t,y,p,tU,ut)
    E_B = y(1); E_L = y(2); A_B = y(3); A_L = y(4);

    % u(t): mol/min input into blood
    u = interp1(tU,ut,t,"previous",0);

    % Units: Q in L/min; V in L; concentrations in mol/L (M)
    Q  = p.QH;
    VL = p.VL;
    VB = p.VB;

    % Reaction rate terms
    vADH  = p.alpha * p.VmaxE * (E_L/(p.KmE + E_L));
    vALDH = p.alpha * p.VmaxA * (A_L/(p.KmA + A_L));

    % ODEs
    dE_B = (Q/VB)*(E_L - E_B) + u/VB;
    dE_L = (Q/VL)*(E_B - E_L) - vADH;

    dA_B = (Q/VB)*(A_L - A_B);
    dA_L = (Q/VL)*(A_B - A_L) + vADH - vALDH;

    dydt = [dE_B; dE_L; dA_B; dA_L];
end
