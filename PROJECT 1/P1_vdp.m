% Pranav Nyati
% 20CS30037
% Proj-1

% Van-der-Pol (VDP) equation:
% y" - u(1- y^2)*y' + y = 0 
% reduction to two first order diff equations:
% let y1 = (1/u)*y', then
% Eqn1 => y1 = (1/u)*y'
% Eqn2 => y1' = u(1-y^2)y1 - y/u

function P1_vdp()
    T_span_1 = [0 100];
    T_span_2 = [0 300];
    y1_0 = 1;
    y2_0 = 0;
    Y_0 = [y1_0; y2_0];
    
    % ODE simulation using ode45 for u = 1
    [T1, Y1] = ode45(@vdp_1, T_span_1, Y_0);

    % ODE simulation using ode45 for u = 0.1
    [T2, Y2] = ode45(@vdp_2, T_span_1, Y_0);

    % ODE simulation using ode45 for u = 100
    tStart_ode45 = tic;
    [T3, Y3] = ode45(@vdp_3, T_span_2, Y_0);
    tEnd_ode45 = toc(tStart_ode45);

    % ODE simulation using ode15 for u = 100 (stiff case)
    tStart_ode15 = tic;
    [T4, Y4] = ode15s(@vdp_3, T_span_2, Y_0);
    tEnd_ode15 = toc(tStart_ode15);

    % Comparing time taken by ODE45 and ODE15 for u = 100 case
    fprintf("Time taken by ODE45 for u = 100: %f milliseconds\n", tEnd_ode45*(10^3));
    fprintf("Time taken by ODE15s for u = 100: %f milliseconds\n", tEnd_ode15*(10^3));

    % y-t plots for different u
    figure()
    subplot(3, 1, 1); plot(T2, Y2(:, 1), T2, Y2(:, 2)); 
    legend('y_1','y_2');title("y-t plot for u = 0.1");
    subplot(3, 1, 2); plot(T1, Y1(:, 1), T1, Y1(:, 2)); 
    legend('y_1','y_2');title("y-t plot for u = 1");
    subplot(3, 1, 3); plot(T3, Y3(:, 1), T3, Y3(:, 2)); 
    legend('y_1','y_2');title("y-t plot for u = 100 (ode45)");

    % phase-plane analysis plots for the 4 cases same as above
    figure()
    subplot(3, 1, 1); plot(Y2(:, 1), Y2(:, 2)); title("phase plane plot for u = 0.1");
    subplot(3, 1, 2); plot(Y1(:, 1), Y1(:, 2)); title("phase plane plot for u = 1");
    subplot(3, 1, 3); plot(Y3(:, 1), Y3(:, 2)); title("phase plane plot for u = 100 (ode45)");
end

% VDP for u = 1
function der_vec = vdp_1(t, y)
    der_vec = [y(2); (1- (y(1)^2))*y(2) - y(1)];
end

% VDP for u = 0.1
function der_vec = vdp_2(t, y)
    der_vec = [0.1*y(2); 0.1*(1- (y(1)^2))*y(2) - (y(1)/(0.1))];
end

% VDP for u = 100
function der_vec = vdp_3(t, y)
    der_vec = [100*y(2); 100*(1- (y(1)^2))*y(2) - (y(1)/(100))];
end

