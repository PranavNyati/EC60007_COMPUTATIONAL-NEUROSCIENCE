% Pranav Nyati
% 20CS30037
% Proj-2


%%
function COMP_NEURO_P2_20CS30037()
    
    fprintf("\n\n  MORIS LECAR EQUATIONS\n\n")
    
    % QUESTION 2 - Plotting null-clines, finding equilibrium points and
    % their eigenvalues, quiver plots
    
    fprintf("\nQUESTION 2:\n\n");
    % First Set of Parameters for MLE:
    g_Ca = 4.4; g_K = 8.0; g_L = 2.0;
    V_Ca = 120; V_K = -84; V_L = -60; phi = 0.02;
    V1 = -1.2; V2 = 18; V3 = 2; V4 = 30; C = 20; I_ext = 0.0;

    V_temp = linspace(-80, 60, 400);
 

    % Method 1 of calculation:
    figure(1)
    hold on
    
    % w-nullcline 
    % dw/dt = 0 => w = w_inf(V)
    w_nullcline = 0.5*(1 + tanh((V_temp - V3)/V4));

    m_inf_V = 0.5*(1 + tanh((V_temp - V1)/V2));

    % v-nullcline
    % dV/dt = 0
    Ca_term = g_Ca*(V_temp - V_Ca);
    K_term = g_K*(V_temp - V_K);
    L_term = g_L*(V_temp - V_L);

    v_nullcline = (I_ext - (m_inf_V.*Ca_term) - L_term)./(K_term);

    plot1 = plot(V_temp, min(100*w_nullcline, 100), 'Color', 'red', 'LineStyle','--', LineWidth=2);
    plot2 = plot(V_temp, min(100*v_nullcline, 100), 'Color', 'blue', 'LineStyle','--', 'LineWidth',2);

    yline(0, 'Color', 'black', 'LineWidth', 1)
    title('Nullclines and Equilibrium points for MLE Parameter Set 1');
    xlabel('V');
    ylabel('100*W');
    
    
    % finding the equilibrium point by iterating through the different
    % values of V

    for i = 1:numel(V_temp)
        p1 = w_nullcline(i);
        p2 = v_nullcline(i);

        if (abs(p2 - p1) < 0.002)
            V_eq = V_temp(i);
            W_eq = (p2 + p1)/2;
            fprintf("Equilibrium point for parameter set 1 for the case of I_ext = 0 for MLE:\n");
            fprintf("W_eq = %.3f, V_eq = %.3f\n", W_eq, V_eq);

        end
    end

  
    hold on
    plot3 = plot(V_eq, 100*W_eq, 'Marker','o', 'Color', 'black', 'MarkerFaceColor','green');
    
    % Method-2: Finding the equilibrium points using sym package of Matlab

    F = @(x) [(1/C)*(I_ext + (-g_Ca*( (0.5*( 1 + tanh((x(1)-V1)/(V2))))*(x(1)-V_Ca))) + (-g_K*(x(2)*(x(1)-V_K) )) + (-g_L*(x(1) - V_L)));  phi*(0.5*(1 + tanh((x(1)-V3)/(V4)) ) - x(2))/(1/cosh((x(1)-V3)/(2*V4)))];
    starting_pt = [-60; 0.01];
    [x,~] = fsolve(F,starting_pt);
    disp("Equilibrium point for MLE with first set of variables, i external = 0")
    Eq_pt_method2 = x;
    disp(Eq_pt_method2);
        
    V_eq = Eq_pt_method2(1); W_eq = Eq_pt_method2(2);

    % plotting the quiver plots

    [v_arrow,w_arrow] = meshgrid(linspace(-80,60,35), linspace(0,1,35));
    m_inf_v_arrow = 0.5*(1 + tanh((v_arrow-V1)/(V2))); 

    tau_w = 1./cosh((v_arrow-V3)/(2*V4));
    dv_dt = (1/C)*(I_ext - g_Ca*(m_inf_v_arrow.*(v_arrow-V_Ca)) - g_K*(w_arrow.*(v_arrow-V_K)) - g_L*(v_arrow - V_L));
    dw_dt = phi * (0.5*( 1 + tanh((v_arrow-V3)/(V4)) ) - w_arrow)./tau_w;
    hold on
    plot4 = quiver(v_arrow,100*w_arrow, dv_dt, 100*dw_dt, 3, 'color',[0 0 0]); % arrow length scaled 2 times for visibility

    legend([plot1, plot2, plot3, plot4], ["W-nullcline", "V-nullcline", "Equilibrium Point", "Gradient Direction"]);    

    grid on
    hold off
    

    % Question 3: Nature of Equilibrium point
    
    fprintf("\nQUESTION 3:\n\n");
    % Evaluating the Jacobian at the Equilibrium point
    
    syms V_temp W_temp
    % f1 = dv/dt; f2 = dw/dt
    f1 = (1/C)*(I_ext + (-g_Ca*((0.5*(1 + tanh((V_temp-V1)/(V2)) ))*(V_temp-V_Ca) )) + (-g_K*(W_temp*(V_temp-V_K))) + (-g_L*(V_temp-V_L)));
    f2 = phi*(0.5*(1 + tanh((V_temp-V3)/(V4)) ) - W_temp)/(1/cosh((V_temp-V3)/(2*V4)));

    %Jacobian Matrix
    j_mat = jacobian([f1; f2], [V_temp; W_temp]);
    j_mat_at_eqlbrm = subs(j_mat, {sym('V_temp'), sym('W_temp')}, {V_eq, W_eq});

    fprintf("\nJacobian matrix for the first set of MLE parameters (for I_ext = 0): \n");
    disp(double(j_mat_at_eqlbrm));
    fprintf("\n");
    eigen_vals = double(eig(j_mat_at_eqlbrm));
    fprintf("The eigenvalues for the Jacobian matrix are:\n");
    disp(eigen_vals);
    fprintf("\n");
    
    res = check_stability(eigen_vals(1), eigen_vals(2));
    fprintf("The given equilibrium point is a %s.\n", res);


    % Question 5: Action potential using different values of phi

    T_span = [0 500]; 
    y_1_0 = V_eq + 46.5;
    y_2_0 = W_eq;
    Y_0 = [y_1_0; y_2_0];     % initial value
    
    % Simulation of MLE with parameter set 1 for phi_1 = 0.02 using ODE15s
    [~, Y1] = ode15s(@MLE_phi_1, T_span, Y_0);

    % Simulation of MLE with parameter set 1 for phi_2 = 0.04 using ODE15s
    [~, Y2] = ode15s(@MLE_phi_2, T_span, Y_0);

    % Simulation of MLE with parameter set 1 for phi_3 = 0.01 using ODE15s
    [~, Y3] = ode15s(@MLE_phi_3, T_span, Y_0);

    
    % Plotting the three action potentials corresponding to different
    % values of phi
    figure(2)
    hold on
    grid on

    % w-nullcline 
    V_temp = linspace(-80, 60, 400);

    w_nullcline = 0.5*(1 + tanh((V_temp - V3)/V4));
    m_inf_V = 0.5*(1 + tanh((V_temp - V1)/V2));

    % v-nullcline
    Ca_term = g_Ca*(V_temp - V_Ca);
    K_term = g_K*(V_temp - V_K);
    L_term = g_L*(V_temp - V_L);
    v_nullcline = (I_ext - (m_inf_V.*Ca_term) - L_term)./(K_term);

    % plotting the nullclines 
    plot1 = plot(V_temp, max(min(100*w_nullcline, 100), 0), 'Color', 'red', 'LineStyle','-', 'LineWidth', 1);
    plot2 = plot(V_temp, max(min(100*v_nullcline, 100), 0), 'Color', 'blue', 'LineStyle','-', 'LineWidth', 1);

    yline(0, 'Color', 'black', 'LineWidth', 1)
    
    % plot the equilibrium point
    hold on
    plot3 = plot(V_eq, 100*W_eq, 'Marker','o', 'Color', 'black', 'MarkerFaceColor', 'green');
    plot8 = plot(y_1_0, 100*y_2_0, 'Marker','o', 'Color', 'black', 'MarkerFaceColor','blue');

    % plotting the phase plane of the Action Potentials for different phi
    plot4 = plot(Y1(:, 1), (100*Y1(:, 2)));
    hold on
    plot5 = plot(Y2(:, 1), (100*Y2(:, 2)));
    hold on
    plot6 = plot(Y3(:, 1), (100*Y3(:, 2)));
    
    hold on
    % quiver plot
    [v_arrow,w_arrow] = meshgrid(linspace(-80,60,20), linspace(0,1,20));
    m_inf_v_arrow = 0.5*(1 + tanh((v_arrow-V1)/(V2))); 

    tau_w = 1./cosh((v_arrow-V3)/(2*V4));
    dv_dt = (1/C)*(I_ext - g_Ca*(m_inf_v_arrow.*(v_arrow-V_Ca)) - g_K*(w_arrow.*(v_arrow-V_K)) - g_L*(v_arrow - V_L));
    dw_dt = phi * (0.5*( 1 + tanh((v_arrow-V3)/(V4)) ) - w_arrow)./tau_w;
    hold on
    plot7 = quiver(v_arrow,100*w_arrow, dv_dt, 100*dw_dt, 2, 'color',[0 0 0]); 

    title('Phase Plane plots of Action Potentials for different phi');
    xlabel('V');
    ylabel('100*W');
    legend([plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8], ["W-nullcline", "V-nullcline", "Equilibrium Point", "AP phi_1=0.02", "AP phi_2=0.04", "AP phi_3 = 0.01", "Gradient Direction", "Starting Point"]);
    hold off
   

    % Question 6: Depolarizing current pulses of different amplitudes 

    figure(3)
    hold on
    grid on

    % w-nullcline 
    V_temp = linspace(-80, 60, 400);

    w_nullcline = 0.5*(1 + tanh((V_temp - V3)/V4));
    m_inf_V = 0.5*(1 + tanh((V_temp - V1)/V2));

    % v-nullcline
    Ca_term = g_Ca*(V_temp - V_Ca);
    K_term = g_K*(V_temp - V_K);
    L_term = g_L*(V_temp - V_L);
    v_nullcline = (I_ext - (m_inf_V.*Ca_term) - L_term)./(K_term);

    % plotting the nullclines 
    plot1 = plot(V_temp, max(min(100*w_nullcline, 100), 0), 'Color', 'red', 'LineStyle','-', 'LineWidth', 1);
    plot2 = plot(V_temp, max(min(100*v_nullcline, 100), 0), 'Color', 'blue', 'LineStyle','-', 'LineWidth', 1);

    yline(0, 'Color', 'black', 'LineWidth', 1)
    
    % plot the equilibrium point
    hold on
    plot3 = plot(V_eq, 100*W_eq, 'Marker','o', 'Color', 'black', 'MarkerFaceColor', 'green');
   

    hold on
     % quiver plot
    [v_arrow,w_arrow] = meshgrid(linspace(-80,60,20), linspace(0,1,20));
    m_inf_v_arrow = 0.5*(1 + tanh((v_arrow-V1)/(V2))); 

    tau_w = 1./cosh((v_arrow-V3)/(2*V4));
    dv_dt = (1/C)*(I_ext - g_Ca*(m_inf_v_arrow.*(v_arrow-V_Ca)) - g_K*(w_arrow.*(v_arrow-V_K)) - g_L*(v_arrow - V_L));
    dw_dt = phi * (0.5*( 1 + tanh((v_arrow-V3)/(V4)) ) - w_arrow)./tau_w;
    hold on
    plot4 = quiver(v_arrow,100*w_arrow, dv_dt, 100*dw_dt, 2, 'color',[0 0 0]);
    hold on

    V_shift = [20 30 40 44 45 46 47 48 49 50 60 70];
    V_amp = zeros(1013, 2);
    count = 1;

    for i = 1:numel(V_shift)
        T_span = [0 500]; 
        y_1_0 = V_eq + V_shift(i);
        y_2_0 = W_eq;
        Y_0 = [y_1_0; y_2_0];

        [~, Y] = ode15s(@MLE_phi_1, T_span, Y_0);
        plot(Y(:, 1), (100*Y(:, 2)));
        V_amp(count, 1) = y_1_0;
        V_amp(count, 2) = max(Y(:, 1));
        count = count + 1;
        hold on

    end

   
    for i = 45:0.001:46
        y_1_0 = V_eq + i;
        y_2_0 = W_eq;
        Y_0 = [y_1_0; y_2_0];

        [~, Y] = ode15s(@MLE_phi_1, T_span, Y_0);
        if (mod(100*i, 2) == 0)
            plot(Y(:, 1), (100*Y(:, 2)));
            hold on
        end

        V_amp(count, 1) = y_1_0;
        V_amp(count, 2) = max(Y(:, 1));
        count = count + 1;

    end

    title('Threshold behaviour of Neurons depicted by MLE');
    xlabel('V');
    ylabel('100*W');
    legend([plot1, plot2, plot3, plot4], ["W-nullcline", "V-nullcline", "Equilibrium Point", "Gradient Direction"]);

    hold off
  
    V_amp = unique(V_amp(:, 1:2), 'rows');

    % plot of maximum voltage amplitude versus initial voltage values to
    % depict the threshold behaviour

    figure(4)
    hold on
    grid on
    plot(V_amp(:, 1), V_amp(:, 2), Marker="o", MarkerEdgeColor="black", MarkerSize=3, MarkerFaceColor="auto");
  
    title('Plot of Maximum Voltage Amplitude for different Initial Voltages (presence of Threshold Behaviour)');
    xlabel('Initial Voltage');
    ylabel('Maximum Voltage Amplitude');
    hold off

    % Question 7:  Steady external current for three different steady
    % states

    fprintf("\nQUESTION 7:\n\n");

    figure(5)
    hold on
    grid on

    I_ext = 86.00;

     % w-nullcline 
    V_temp = linspace(-74, 60, 400);

    w_nullcline = 0.5*(1 + tanh((V_temp - V3)/V4));
    m_inf_V = 0.5*(1 + tanh((V_temp - V1)/V2));

    % v-nullcline
    Ca_term = g_Ca*(V_temp - V_Ca);
    K_term = g_K*(V_temp - V_K);
    L_term = g_L*(V_temp - V_L);
    v_nullcline = (I_ext - (m_inf_V.*Ca_term) - L_term)./(K_term);

    % plotting the nullclines 
    plot1 = plot(V_temp, max(min(100*w_nullcline, 200), 0), 'Color', 'red', 'LineStyle','--', 'LineWidth', 1.5);
    plot2 = plot(V_temp, max(min(100*v_nullcline, 200), 0), 'Color', 'blue', 'LineStyle','--', 'LineWidth', 1.5);

    % Plotting the y = 0 line
    yline(0, 'Color', 'black', 'LineWidth', 1)


    % finding the equilibrium point by iterating through the different
    % values of V for the steady state current case in MLE

    for i = 1:numel(V_temp)
        p1 = w_nullcline(i);
        p2 = v_nullcline(i);

        if (abs(p2 - p1) < 0.001)
            V_eq_steady_i = V_temp(i);
            W_eq_steady_i = (p2 + p1)/2;
            fprintf("\nEquilibrium point for the case of Steady Externel Current Density = 86 uA/cm2 for MLE:\n");
            fprintf("W_eq = %.3f, V_eq = %.3f\n", W_eq_steady_i, V_eq_steady_i);

        end
    end

    % plot the equilibrium point of the steady current case
    hold on
    plot3 = plot(V_eq_steady_i, 100*W_eq_steady_i, 'Marker','o', 'Color', 'black', 'MarkerFaceColor', 'magenta');

    % Plotting the quiver plot
    hold on

    [v_arrow,w_arrow] = meshgrid(linspace(-70,60,30), linspace(0,1.5,30));
    m_inf_v_arrow = 0.5*(1 + tanh((v_arrow-V1)/(V2))); 

    tau_w = 1./cosh((v_arrow-V3)/(2*V4));
    dv_dt = (1/C)*(I_ext - g_Ca*(m_inf_v_arrow.*(v_arrow-V_Ca)) - g_K*(w_arrow.*(v_arrow-V_K)) - g_L*(v_arrow - V_L));
    dw_dt = phi * (0.5*( 1 + tanh((v_arrow-V3)/(V4)) ) - w_arrow)./tau_w;
    plot4 = quiver(v_arrow,100*w_arrow, dv_dt, 100*dw_dt, 2, 'color',[0 0 0]);

    % Plotting the three trajectories with 3 different initial conditions

    T_span = [0 1000]; 

    % Initial condition 1: Equilibrium point of original case (no steady
    % current)
    
    y_i1_1 = V_eq;
    y_i1_2 = W_eq;
    Y_init_1 = [y_i1_1; y_i1_2];     % initial value

    [~, Y1] = ode15s(@MLE_steady_current, T_span, Y_init_1);

    % Initial condition 2: Equilibrium point of the steady current case 
    
    y_i2_1 = V_eq;
    y_i2_2 = W_eq;
    Y_init_2 = [y_i2_1; y_i2_2];     % initial value

    % Simulation of MLE with steady state current and initial condition 1
    [~, Y2] = ode15s(@MLE_steady_current, T_span, Y_init_2);

    % Initial condition 3: V = -27.9, W = 0.17
    
    y_i3_1 = -27.9;
    y_i3_2 = 0.17;
    Y_init_3 = [y_i3_1; y_i3_2];     % initial value

     % Simulation of MLE with steady state current and initial condition 1
    [~, Y3] = ode15s(@MLE_steady_current, T_span, Y_init_3);

    hold on
    plot5 = plot(Y1(:, 1), (100*Y1(:, 2)), LineWidth=1);
    hold on
    plot6 = plot(Y2(:, 1), (100*Y2(:, 2)), 'Color', 'green', LineWidth=1);
    hold on
    plot7 = plot(Y3(:, 1), (100*Y3(:, 2)), LineWidth=1);
    
    % Marking the 3 different starting points/initial conditions
    plot8 = plot(V_eq, 100*W_eq, 'Marker','o', 'Color', 'black', 'MarkerFaceColor', 'green');
    plot9 = plot(y_i3_1, 100*y_i3_2, 'Marker','o', 'Color', 'black', 'MarkerFaceColor', 'blue');


    title('Phase Plane Plots and Trajectories for MLE with Steady State Current for 3 different Initial Conditions');
    xlabel('V');
    ylabel('100*W');
    legend([plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9], ["W-nullcline", "V-nullcline","Equilibrium Point & Initial Condition 2", "Gradient Direction", "Trajectory 1", "Trajectory 2", "Trajectory 3", "Initial Condition 1", "Initial Condition 3"]);
    hold off
  
    % Question 8: To find the contour dividing the phase plane into the two
    % different stable states (limit cycle and equilibrium pt)
    % To find an unstable periodic orbit like UPO by running the system
    % backwards in time

    figure(6)
    hold on
    grid on

    % plotting the nullclines 
    plot1_fig6 = plot(V_temp, max(min(100*w_nullcline, 200), 0), 'Color', 'red', 'LineStyle','--', 'LineWidth', 1.5);
    plot2_fig6 = plot(V_temp, max(min(100*v_nullcline, 200), 0), 'Color', 'blue', 'LineStyle','--', 'LineWidth', 1.5);

    % Plotting the y = 0 line
    yline(0, 'Color', 'black', 'LineWidth', 1)
    
    % plot the equilibrium point of the steady current case
    hold on
    plot3_fig6 = plot(V_eq_steady_i, 100*W_eq_steady_i, 'Marker','o', 'Color', 'black', 'MarkerFaceColor', 'magenta');

    % Plotting the quiver plot
    hold on
    plot4_fig6 = quiver(v_arrow,100*w_arrow, dv_dt, 100*dw_dt, 2, 'color',[0 0 0]);

    % Plotting the UPO (contour separating the points converging to stable
    % limit cycle with the points converging to the stable equilibrium at
    % the intersection of nullclines

    T_span = [0 1000]; 

    % On reversing time, the stable equilibrium point changes to unstable,
    % so starting from slightly away from it, we can simulate the system to
    % obtain the UPO
    y_i1_1 = V_eq_steady_i + 0.1;
    y_i1_2 = W_eq_steady_i + 0.001;
    Y_init_1 = [y_i1_1; y_i1_2];     % initial value

    % Simulation of MLE with steady current and time reversal to get the
    % trajectory of UPO

    [~, Y1] = ode15s(@MLE_steady_current_time_reverse, T_span, Y_init_1);
    
    hold on
    plot5_fig6 = plot(Y1(:, 1), (100*Y1(:, 2)), LineWidth=2, LineStyle="-", Color=[0 0 0]);

    % Another initial condition inside the UPO (so it converges to the
    % stable equilibrium point)
    y_i2_1 = -27.9;
    y_i2_2 = 0.17;
    Y_init_2 = [y_i2_1; y_i2_2];     % initial value

    [~, Y2] = ode15s(@MLE_steady_current, T_span, Y_init_2);

    hold on
    plot6_fig6 = plot(Y2(:, 1), (100*Y2(:, 2)), LineWidth=1);

    % Initial condition 3 (equilibrium point of original MLE with no steady current:
    % Outside the UPO (so it converges to the stable limit cycle)

    y_i3_1 = V_eq;
    y_i3_2 = W_eq;
    Y_init_3 = [y_i3_1; y_i3_2];     % initial value

    % Simulation of MLE with steady state current and initial condition 1
    [~, Y3] = ode15s(@MLE_steady_current, T_span, Y_init_3);

    hold on
    plot7_fig6 = plot(Y3(:, 1), (100*Y3(:, 2)), LineWidth=1);

    % plotting the starting points of the above 3 trajectories
    plot8_fig6 = plot(y_i1_1, 100*y_i1_2, 'Marker','o', 'Color', 'black', 'MarkerFaceColor', 'red');
    plot9_fig6 = plot(y_i2_1, 100*y_i2_2, 'Marker','o', 'Color', 'black', 'MarkerFaceColor', 'blue');
    plot10_fig6 = plot(y_i3_1, 100*y_i3_2, 'Marker','o', 'Color', 'black', 'MarkerFaceColor', 'green');
    
    % Plotting more trajectories with start state inside the UPO and hence
    % converge to the stable equilibrium point

    V_offset = [-5.0, 7.5, 11, 2];
    W_offset = [0.05, 0.07, 0.0865, 0.075];

    for i = 1:numel(V_offset)
        Y_init_temp = [V_eq_steady_i + V_offset(i), W_eq_steady_i + W_offset(i)];
        
        [~, Y] = ode15s(@MLE_steady_current, T_span, Y_init_temp);

        hold on
        plot(Y(:, 1), (100*Y(:, 2)), LineWidth=1);
    end

    % Plotting more trajectories with start state outside the UPO and hence
    % converge to stable limit cycle

    V_offset_1 = [-10.0, 35, 11, 40, 30, 3, 32];
    W_offset_1 = [0.05, 0.2, 0.4, -0.1, 0.05, 0.2, 0.55];

    for i = 1:numel(V_offset_1)
        Y_init_temp = [V_eq_steady_i + V_offset_1(i), W_eq_steady_i + W_offset_1(i)];
        
        [~, Y] = ode15s(@MLE_steady_current, T_span, Y_init_temp);

        hold on
        plot(Y(:, 1), (100*Y(:, 2)), LineWidth=1);
    end


    % Plot title and labels
    xlim([-75 40]);
    ylim([0 80]);
    title('Plot depicting Unstable Periodic Orbit for Steady Current = 86 uA/cm2');
    xlabel('V');
    ylabel('100*W');
    legend([plot1_fig6, plot2_fig6, plot3_fig6, plot4_fig6, plot5_fig6, plot6_fig6, plot7_fig6, plot8_fig6, plot9_fig6, plot10_fig6], ["W-nullcline", "V-nullcline","Equilibrium Point", "Gradient Direction", "Trajectory 1 - UPO", "Trajectory 2 - Stable equilibrium convergence", "Trajectory 3 - Limit Cycle Convergence", "Start Pt. of Trajectory 1", "Start Pt. of Trajectory 2", "Start Pt. of Trajectory 3"]);

    hold off

    % Question 9:
    fprintf("\nQUESTION 9:\n\n");

    I_ext_vals = [80.0, 86.0, 90.0];
    eqbm_pts = zeros(3, 2);
    figure(7)
       

    for i = 1:numel(I_ext_vals)

        % w-nullcline 
        V_temp = linspace(-70, 60, 400);
    
        w_nullcline = 0.5*(1 + tanh((V_temp - V3)/V4));
        m_inf_V = 0.5*(1 + tanh((V_temp - V1)/V2));
    
        % v-nullcline
        Ca_term = g_Ca*(V_temp - V_Ca);
        K_term = g_K*(V_temp - V_K);
        L_term = g_L*(V_temp - V_L);
        v_nullcline = (I_ext_vals(i) - (m_inf_V.*Ca_term) - L_term)./(K_term);
    
        % plotting the nullclines 
        subplot(2, 2, i);
        hold on
        grid on
        plot1 = plot(V_temp, max(min(100*w_nullcline, 100), 0), 'Color', 'red', 'LineStyle','-', 'LineWidth', 1);
        hold on
        plot2 = plot(V_temp, max(min(100*v_nullcline, 100), 0), 'Color', 'blue', 'LineStyle','-', 'LineWidth', 1);
        hold on
        yline(0, 'Color', 'black', 'LineWidth', 1)
        
        
        % finding the equilibrium point by iterating through the different
        % values of V
    
        for j = 1:numel(V_temp)
            p1 = w_nullcline(j);
            p2 = v_nullcline(j);
    
            if (abs(p2 - p1) < 0.001)
                eqbm_pts(i, 1) = V_temp(j);
                eqbm_pts(i, 2) = (p2 + p1)/2;
                fprintf("\nEquilibrium point for parameter set 1 for the case of I_ext = %f uA/cm2 for MLE:\n", I_ext_vals(i));
                fprintf("W_eq = %.3f, V_eq = %.3f\n\n", eqbm_pts(i, 2), eqbm_pts(i, 1));
    
            end
        end
        
        % plot the equilibrium point
        hold on
        plot3 = plot(eqbm_pts(i, 1), 100*eqbm_pts(i, 2), 'Marker','o', 'Color', 'black', 'MarkerFaceColor', 'green');
        
       
        % Plotting the quiver plot
        hold on
        [v_arrow,w_arrow] = meshgrid(linspace(-70,60,20), linspace(0,1.5,20));
        m_inf_v_arrow = 0.5*(1 + tanh((v_arrow-V1)/(V2))); 
    
        tau_w = 1./cosh((v_arrow-V3)/(2*V4));
        dv_dt = (1/C)*(I_ext_vals(i) - g_Ca*(m_inf_v_arrow.*(v_arrow-V_Ca)) - g_K*(w_arrow.*(v_arrow-V_K)) - g_L*(v_arrow - V_L));
        dw_dt = phi * (0.5*( 1 + tanh((v_arrow-V3)/(V4)) ) - w_arrow)./tau_w;
        plot4 = quiver(v_arrow,100*w_arrow, dv_dt, 100*dw_dt, 2, 'color',[0 0 0]);

        % Labelling the plot
        str_temp = sprintf('Equilibrium points for I external = %f uA/cm^2', I_ext_vals(i));
        title(str_temp);
        legend([plot1, plot2, plot3, plot4], ["W-nullcline", "V-nullcline","Equilibrium Point", "Gradient Direction"]);
        xlabel('V');
        ylabel('100*W');

        hold off

        % Stability Analysis of the Equilibrium Point
        
        syms V_temp W_temp
        % f1 = dv/dt; f2 = dw/dt
        f1 = (1/C)*(I_ext_vals(i) + (-g_Ca*((0.5*(1 + tanh((V_temp-V1)/(V2)) ))*(V_temp-V_Ca) )) + (-g_K*(W_temp*(V_temp-V_K))) + (-g_L*(V_temp-V_L)));
        f2 = phi*(0.5*(1 + tanh((V_temp-V3)/(V4)) ) - W_temp)/(1/cosh((V_temp-V3)/(2*V4)));
    
        %Jacobian Matrix
        j_mat = jacobian([f1; f2], [V_temp; W_temp]);
        j_mat_at_eqlbrm = subs(j_mat, {sym('V_temp'), sym('W_temp')}, {eqbm_pts(i, 1), eqbm_pts(i, 2)});
    
        fprintf("\nJacobian matrix for the first set of MLE parameters (for I_ext = %f uA/cm2): \n", I_ext_vals(i));
        disp(double(j_mat_at_eqlbrm));
        eigen_vals = double(eig(j_mat_at_eqlbrm));
        fprintf("The eigenvalues for the Jacobian matrix are: \n");
        disp(eigen_vals);
        
        res = check_stability(eigen_vals(1), eigen_vals(2));
        fprintf("The given equilibrium point is a %s.\n", res);

    end

    % Rate of firing of action potential vs applied current plot for param
    % set 1

    figure(8)
    hold on
    grid on

    I_ext_arr = 80:0.03:100;
    AP_firing_rates = zeros(numel(I_ext_arr), 1);

    for i = 1:numel(I_ext_arr)

        % Finding Equilibrium point for each value of I_ext
        F = @(x) [(1/C)*(I_ext_arr(i) + (-g_Ca*( (0.5*( 1 + tanh((x(1)-V1)/(V2))))*(x(1)-V_Ca))) + (-g_K*(x(2)*(x(1)-V_K) )) + (-g_L*(x(1) - V_L)));  phi*(0.5*(1 + tanh((x(1)-V3)/(V4)) ) - x(2))/(1/cosh((x(1)-V3)/(2*V4)))];
        starting_pt = [-65; 0.01];
        options = optimset('Display','off');
        [x,~] = fsolve(F,starting_pt, options);

        V_eq_temp = x(1); W_eq_temp = x(2);
        [T, Y] = MLE_for_given_steady_I_param_set1(I_ext_arr(i), V_eq_temp + 0.2, W_eq_temp+ 0.005);

        AP_firing_rates(i) = calculate_firing_rate(T, Y);

    end
  
    % plot for firing rates vs I_ext for range 80-100 uA/cm2
    plot(I_ext_arr, AP_firing_rates, LineWidth=2);
    title('Rate of AP versus applied I ext in range 80-100 uA/cm2 for MLE parameter set 1')
    xlabel('I external (in uA/cm2')
    ylabel('Frequency of Action Potentials per second')
    hold off

    % Question 10: MLE with different set of parameters

    fprintf("\nQUESTION 10:\n\n");
    
    g_Ca = 4.0; g_K = 8.0; g_L = 2.0;
    V_Ca = 120; V_K = -84; V_L = -60; phi = 0.0667;
    V1 = -1.2; V2 = 18; V3 = 12; V4 = 17.4; C = 20; I_ext = 30.0;

    figure(9)
    hold on
    grid on

    % w-nullcline 
    V_temp = linspace(-75, 50, 400);

    w_nullcline = 0.5*(1 + tanh((V_temp - V3)/V4));
    m_inf_V = 0.5*(1 + tanh((V_temp - V1)/V2));

    % v-nullcline
    Ca_term = g_Ca*(V_temp - V_Ca);
    K_term = g_K*(V_temp - V_K);
    L_term = g_L*(V_temp - V_L);
    v_nullcline = (I_ext - (m_inf_V.*Ca_term) - L_term)./(K_term);

    % plotting the nullclines 
    plot1 = plot(V_temp, max(min(100*w_nullcline, 100), 0), 'Color', 'red', 'LineStyle','-', 'LineWidth', 1);
    plot2 = plot(V_temp, max(min(100*v_nullcline, 100), 0), 'Color', 'blue', 'LineStyle','-', 'LineWidth', 1);

    hold on
    yline(0, 'Color', 'black', 'LineWidth', 1);

    % quiver plot
    [v_arrow,w_arrow] = meshgrid(linspace(-70,40,20), linspace(0,1,20));
    m_inf_v_arrow = 0.5*(1 + tanh((v_arrow-V1)/(V2))); 

    tau_w = 1./cosh((v_arrow-V3)/(2*V4));
    dv_dt = (1/C)*(I_ext - g_Ca*(m_inf_v_arrow.*(v_arrow-V_Ca)) - g_K*(w_arrow.*(v_arrow-V_K)) - g_L*(v_arrow - V_L));
    dw_dt = phi * (0.5*( 1 + tanh((v_arrow-V3)/(V4)) ) - w_arrow)./tau_w;
    hold on
    plot3 = quiver(v_arrow,100*w_arrow, dv_dt, 100*dw_dt, 2, 'color',[0 0 0]);
   
    % Finding the equilibrium points:
    F = @(x) [(1/C)*(I_ext + (-g_Ca * ( (0.5 * ( 1 + tanh((x(1)-V1)/(V2)) ))*(x(1)-V_Ca) )) + (-g_K*( x(2)*(x(1)-V_K) )) + (-g_L * (x(1) - V_L)));  phi * (0.5 * ( 1 + tanh((x(1)-V3)/(V4)) ) - x(2))/(1/cosh((x(1)-V3)/(2*V4)))];
    options = optimset('Display','off');

    % Equilibrium point 1
    starting_pt = [-50; 0.01];
    [x,~] = fsolve(F,starting_pt, options);
    eq_pt1 = x;
    fprintf("The first equilibrium point is: \n");
    disp(eq_pt1);
  
    % Equilibrium point 2
    starting_pt = [-30; 0.02];
    [x,~] = fsolve(F,starting_pt, options);
    eq_pt2 = x;
    fprintf("The second equilibrium point is: \n");
    disp(eq_pt2);
    
    % Equilibrium point 3
    starting_pt = [5; 0.1];
    [x,~] = fsolve(F,starting_pt, options);
    eq_pt3 = x;
    fprintf("The third equilibrium point is: \n");
    disp(eq_pt3);

    % Plotting the equilibrium points on the plot
    plot4 = plot(eq_pt1(1), 100*eq_pt1(2), 'Marker','o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red');
    plot5 = plot(eq_pt2(1), 100*eq_pt2(2), 'Marker','o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'magenta');
    plot6 = plot(eq_pt3(1), 100*eq_pt3(2), 'Marker','o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green');
 

    % Stability analysis of the equilibrium points:

    syms V_temp W_temp
    % f1 = dv/dt; f2 = dw/dt
    f1 = (1/C)*(I_ext + (-g_Ca*((0.5*(1 + tanh((V_temp-V1)/(V2)) ))*(V_temp-V_Ca) )) + (-g_K*(W_temp*(V_temp-V_K))) + (-g_L*(V_temp-V_L)));
    f2 = phi*(0.5*(1 + tanh((V_temp-V3)/(V4)) ) - W_temp)/(1/cosh((V_temp-V3)/(2*V4)));

    %Jacobian Matrix
    j_mat = jacobian([f1; f2], [V_temp; W_temp]);

    % Equilibrium point 1:
    j_mat_at_eqlbrm_1 = subs(j_mat, {sym('V_temp'), sym('W_temp')}, {eq_pt1(1), eq_pt1(2)});

    fprintf("\nJacobian matrix for Equilibrium Point 1 for the 2nd set of MLE parameters (for I_ext = %f uA/cm2): \n", I_ext);
    disp(double(j_mat_at_eqlbrm_1));
    eigen_vals_1 = double(eig(j_mat_at_eqlbrm_1));
    fprintf("The eigenvalues for the Jacobian matrix are: \n");
    disp(eigen_vals_1);

    res1 = check_stability(eigen_vals_1(1), eigen_vals_1(2));
    fprintf("The given equilibrium point is a %s.\n", res1);

    % Equilibrium point 2:
    j_mat_at_eqlbrm_2 = subs(j_mat, {sym('V_temp'), sym('W_temp')}, {eq_pt2(1), eq_pt2(2)});

    fprintf("\nJacobian matrix for Equilibrium Point 2 for the 2nd set of MLE parameters (for I_ext = %f uA/cm2): \n", I_ext);
    disp(double(j_mat_at_eqlbrm_2));
    [eigen_vals_2]  = double(eig(j_mat_at_eqlbrm_2));
    fprintf("The eigenvalues for the Jacobian matrix are: \n");
    disp(eigen_vals_2);
    
    res2 = check_stability(eigen_vals_2(1), eigen_vals_2(2));
    fprintf("The given equilibrium point is a %s.\n", res2);

    % Equilibrium point 3:
    j_mat_at_eqlbrm_3 = subs(j_mat, {sym('V_temp'), sym('W_temp')}, {eq_pt3(1), eq_pt3(2)});

    fprintf("\nJacobian matrix for Equilibrium Point 3 for the 2nd set of MLE parameters (for I_ext = %f uA/cm2): \n", I_ext);
    disp(double(j_mat_at_eqlbrm_3));
    eigen_vals_3 = double(eig(j_mat_at_eqlbrm_3));
    fprintf("The eigenvalues for the Jacobian matrix are: \n");
    disp(eigen_vals_3);
    
    res3= check_stability(eigen_vals_3(1), eigen_vals_3(2));
    fprintf("The given equilibrium point is a %s.\n", res3);


    % Plotting the stable and unstable manifolds of the saddle point
    % (Eq_pt_2)

    [V, ~] = eig(j_mat_at_eqlbrm_2);
 
    neg_eigenval_eigenvec = double(V(:, 1));
    disp(neg_eigenval_eigenvec);
    pos_eigenval_eigenvec = double(V(:, 2));
    disp(pos_eigenval_eigenvec);
    % UNSTABLE MANIFOLDS
    start_pt_unstable_manifold_1 = eq_pt2 + 0.1*(pos_eigenval_eigenvec/(norm(pos_eigenval_eigenvec)));
    start_pt_unstable_manifold_2 = eq_pt2 - 0.1*(pos_eigenval_eigenvec/(norm(pos_eigenval_eigenvec)));

    % Simulating the trajectories of the two unstable manifolds
    [~, Y1] = MLE_for_given_steady_I_param_set2(I_ext, start_pt_unstable_manifold_1(1), start_pt_unstable_manifold_1(2), [0 200]);
    [~, Y2] = MLE_for_given_steady_I_param_set2(I_ext, start_pt_unstable_manifold_2(1), start_pt_unstable_manifold_2(2), [0 200]);



    % Unstable manifold 1
    plot7 = plot(Y1(:, 1), (100*Y1(:, 2)), 'Color', 'green',LineWidth=1.5);
    hold on
    % Unstable manifold 2
    plot8 = plot(Y2(:, 1), (100*Y2(:, 2)), 'Color', 'cyan', LineWidth=1.5);

    % STABLE MANIFOLDS
    
    % Simulating the trajectories of the two stable manifolds by time
    % reversal

    start_pt_unstable_manifold_1 = eq_pt2 - 0.008*(neg_eigenval_eigenvec/(norm(neg_eigenval_eigenvec)));
    %  start_pt_unstable_manifold_1 = eq_pt2 - [0.001,-0.00001];
    start_pt_unstable_manifold_2 = eq_pt2 + 0.01*(neg_eigenval_eigenvec/(norm(neg_eigenval_eigenvec)));

    [~, Y1] = MLE_for_given_steady_I_param_set2(I_ext, start_pt_unstable_manifold_1(1), start_pt_unstable_manifold_1(2), [0 -90]);
    [~, Y2] = MLE_for_given_steady_I_param_set2(I_ext, start_pt_unstable_manifold_2(1), start_pt_unstable_manifold_2(2), [0 -300]);

    
    % Stable manifold 1
    plot9 = plot(Y1(:, 1), (100*Y1(:, 2)), 'Color', 'black', LineWidth=3);
    hold on
    % Stable manifold 2
    plot10 = plot(Y2(:, 1), (100*Y2(:, 2)), LineWidth=1.5);
    
    % Labelling the plot
    str_temp = sprintf('Phase Plane Plot of MLE for 2nd set of Parameters for I ext = %f uA/cm2', I_ext);
    title(str_temp);
    legend([plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10], ["W-nullcline", "V-nullcline", "Gradient Direction", "Equilibrium Point 1 - Stable Point", "Equilibrium Point 2 - Saddle Point", "Equilibrium Point 3 - Unstable Spiral", "Unstable Manifold 1", "Unstable Manifold 2", "Stable Mainfold 1", "Stable Manifold 2"]);
    xlabel('V');
    ylabel('100*W');
    hold off

    %%

    % Question 11:  Change in equilibrium points with varying current
    fprintf("\nQUESTION 11:\n");

    I_arr = 30:38;
    I_arr1 = 39:0.1:40;
    I_arr2 = 41:50;
    I_arr1 = cat(2, I_arr1, I_arr2);
    I_arr = cat(2, I_arr, I_arr1);

    for i = 1:numel(I_arr)
    
        I_ext = I_arr(i);

        % Finding the equilibrium points and stability analysis:

        fprintf("\n\nEQUILIBRIUM POINTS OF I_EXT = %f uA/cm2 => \n", I_ext);
            
        syms V_temp W_temp
        % f1 = dv/dt; f2 = dw/dt
        f1 = (1/C)*(I_ext + (-g_Ca*((0.5*(1 + tanh((V_temp-V1)/(V2)) ))*(V_temp-V_Ca) )) + (-g_K*(W_temp*(V_temp-V_K))) + (-g_L*(V_temp-V_L)));
        f2 = phi*(0.5*(1 + tanh((V_temp-V3)/(V4)) ) - W_temp)/(1/cosh((V_temp-V3)/(2*V4)));
        
        j_mat = jacobian([f1; f2], [V_temp; W_temp]);

        F = @(x) [(1/C)*(I_ext + (-g_Ca * ( (0.5 * ( 1 + tanh((x(1)-V1)/(V2)) ))*(x(1)-V_Ca) )) + (-g_K*( x(2)*(x(1)-V_K) )) + (-g_L * (x(1) - V_L)));  phi * (0.5 * ( 1 + tanh((x(1)-V3)/(V4)) ) - x(2))/(1/cosh((x(1)-V3)/(2*V4)))];
        options = optimset('Display','off');

        % Equilibrium point 1
        starting_pt = [-50; 0.01];
        [x,~] = fsolve(F,starting_pt, options);
        eq_pt1 = x;

        j_mat_at_eqlbrm_1 = subs(j_mat, {sym('V_temp'), sym('W_temp')}, {eq_pt1(1), eq_pt1(2)});
        eigen_vals_1 = double(eig(j_mat_at_eqlbrm_1));
        res1 = check_stability(eigen_vals_1(1), eigen_vals_1(2));
       
        fprintf("The first equilibrium point V_eq = %f, W_eq = %f is a %s.\n", eq_pt1(1), eq_pt1(2), res1);
        
        % Equilibrium point 2
        starting_pt = [-20; 1];
        [x,~] = fsolve(F,starting_pt, options);
        eq_pt2 = x;

        j_mat_at_eqlbrm_2 = subs(j_mat, {sym('V_temp'), sym('W_temp')}, {eq_pt2(1), eq_pt2(2)});
        eigen_vals_2 = double(eig(j_mat_at_eqlbrm_2));
        res2 = check_stability(eigen_vals_2(1), eigen_vals_2(2));

        fprintf("The second equilibrium point V_eq = %f, W_eq = %f is a %s.\n", eq_pt2(1), eq_pt2(2), res2);

        % Equilibrium point 3
        starting_pt = [5; 2];
        [x,~] = fsolve(F,starting_pt, options);
        eq_pt3 = x;
        
        j_mat_at_eqlbrm_3 = subs(j_mat, {sym('V_temp'), sym('W_temp')}, {eq_pt3(1), eq_pt3(2)});
        eigen_vals_3 = double(eig(j_mat_at_eqlbrm_3));
        res3 = check_stability(eigen_vals_3(1), eigen_vals_3(2));

        fprintf("The third equilibrium point V_eq = %f, W_eq = %f is a %s.\n", eq_pt3(1), eq_pt3(2), res3);

    end

    % plot of the rate of firing action potentials versus the applied
    % current for param set 2

    figure(10)
    hold on
    grid on

    I_ext_arr = 30:0.03:45;
    AP_firing_rates = zeros(numel(I_ext_arr), 1);

    for i = 1:numel(I_ext_arr)

        % Finding Equilibrium point for each value of I_ext
        F = @(x) [(1/C)*(I_ext_arr(i) + (-g_Ca*( (0.5*( 1 + tanh((x(1)-V1)/(V2))))*(x(1)-V_Ca))) + (-g_K*(x(2)*(x(1)-V_K) )) + (-g_L*(x(1) - V_L)));  phi*(0.5*(1 + tanh((x(1)-V3)/(V4)) ) - x(2))/(1/cosh((x(1)-V3)/(2*V4)))];
        starting_pt = [-65; 0.01];
        options = optimset('Display','off');
        [x,~] = fsolve(F,starting_pt, options);

        V_eq_temp = x(1); W_eq_temp = x(2);
        [T, Y] = MLE_for_given_steady_I_param_set2(I_ext_arr(i), V_eq_temp + 0.2, W_eq_temp+ 0.005, [0 2000]);

        AP_firing_rates(i) = calculate_firing_rate(T, Y);

    end
  
    % plot for firing rates vs I_ext for range 80-100 uA/cm2
    plot(I_ext_arr, AP_firing_rates, LineWidth=2);
    title('Rate of AP versus applied I ext in range 30-45 uA/cm2 for MLE parameter set 2')
    xlabel('I external (in uA/cm2')
    ylabel('Frequency of Action Potentials per second')
    hold off

    % ###################### HODGKIN HUXLEY EQUATIONS ##################
    
    fprintf("\n\n   HODGKIN HUXLEY EQUATIONS    \n\n");

    % Question 12 and 13: Determination of E_Leak for V_rest = -60 mV
    % Steady state, so dV/dt = 0; dm/dt = dn/dt = dh/dt = 0;

    fprintf("\nQUESTION 12 and 13:\n\n");

    syms V_mem

    G_K_const = 36; G_Na_const = 120; G_L_const = 0.3;
    E_K = -72; E_Na = 55; C = 1;

    alpha_n = -0.01*(V_mem + 50)/(exp(-(V_mem + 50)/10) - 1);
    beta_n = 0.125*exp(-(V_mem + 60)/80);

    alpha_m = -0.1*(V_mem + 35)/(exp(-(V_mem + 35)/10) - 1);
    beta_m = 4*(exp(-(V_mem + 60)/18));

    alpha_h = 0.07*exp(-(V_mem + 60)/20);
    beta_h = 1/(exp(-(V_mem + 30)/10) + 1);

    m_inf_V_rest = alpha_m/(alpha_m + beta_m);
    n_inf_V_rest = alpha_n/(alpha_n + beta_n);
    h_inf_V_rest = alpha_h/(alpha_h + beta_h);

    I_ext = 0.0;

    E_L = V_mem - ((1/G_L_const)*(I_ext + (-G_K_const*(n_inf_V_rest^4)*(V_mem - E_K)) + (-G_Na_const*(m_inf_V_rest^3)*(h_inf_V_rest)*(V_mem - E_Na))));
    E_L_rest = double(subs(E_L, {sym('V_mem')}, {-60}));
    
    fprintf("\nThe value of E_L for V_rest = -60 mV is %f mV\n", E_L_rest);


    % Action Potenitals using HH equations for I_ext = 10 uA/cm2
    
    I_ext_HH = 10;
    T_span = [0 80];
    V_init = -60;
    n_init = double(subs(n_inf_V_rest, {sym('V_mem')}, {V_init}));
    m_init = double(subs(m_inf_V_rest, {sym('V_mem')}, {V_init}));
    h_init = double(subs(h_inf_V_rest, {sym('V_mem')}, {V_init}));


    init_condn = [V_init, n_init, m_init, h_init];

    [T, Y] = HH_equations_for_given_current(I_ext_HH, T_span, init_condn);

    % Action Potential plot
    figure(11)
    
    subplot(2, 1, 1);
    hold on
    grid on
    plot(T, Y(:, 1));

    title('Action Potentials using HH model for I ext = 10 uA/cm2');
    xlabel('Time (in ms)');     
    ylabel('Membrane Potential (in mV)');
    ylim([-80 60]);

    hold off

    % Plot showing variation of n, m, h with time corresponding to above
    % Action Potentials
  

    subplot(2, 1, 2);
    hold on
    grid on

    plot1 = plot(T, Y(:, 2));  
    hold on
    plot2 = plot(T, Y(:, 3));  
    hold on;
    plot3 = plot(T, Y(:, 4));

    title('Variation of n, m, h using HH model for I ext = 10 uA/cm2');
    xlabel('Time (in ms)');     
    ylabel('Ion-Channel Opening Probabilities');
    ylim([0 1]);
    legend([plot1, plot2, plot3], ["n(t)", "m(t)", "h(t)"]);
    hold off
   
    % Question 14:  Model stabiity at Equilibrium for I_ext = 0 at V_mem = -60 mV
    fprintf("\nQUESTION 14:\n\n");

    I_ext = 0;
    syms n m h

    f1 = (1/C)*(I_ext + (-G_K_const*(n^4)*(V_mem - E_K)) + (-G_Na_const*(m^3)*h*(V_mem - E_Na)) + (-G_L_const*(V_mem - E_L_rest)));         % dV/dt
    f2 = alpha_n*(1 - n) - beta_n*n;                                                                                                          % dn/dt
    f3 = alpha_m*(1 - m) - beta_m*m;                                                                                                          % dm/dt
    f4 = alpha_h*(1 - h) - beta_h*h;                                                                                                          % dh/dt
    eqbm = vpasolve([f1 == 0, f2 == 0, f3 == 0, f4 == 0], [V_mem, n, m, h]);

    fprintf("The Equilibrium condition for HH model for I_ext = 0 uA/cm2  and V_rest = -60 mV is:\n");
    fprintf("V_eq = %f mV, n_eq = %f, m_eq = %f, h_eq = %f\n", eqbm.V_mem, eqbm.n, eqbm.m, eqbm.h);
    
    % Checking the stability at equilibrium by starrting slighly away from
    % the Equilibrium position
    init_condn = [double(eqbm.V_mem) + 0.1, double(eqbm.n) + 0.01, double(eqbm.m) + 0.001, double(eqbm.h) + 0.01];
    [T, Y] = HH_equations_for_given_current(I_ext, T_span, init_condn);

    % Action Potential plot
    figure(12)

    subplot(2, 1, 1);
    hold on
    grid on
    plot(T, Y(:, 1));

    title('Action Potentials using HH model for I ext = 0 uA/cm2 starting near Eqbm position');
    xlabel('Time (in ms)');     
    ylabel('Membrane Potential (in mV)');
    ylim([-70 -50]);

    hold off

    % Plot showing variation of n, m, h with time corresponding to above
    % Action Potentials
    
    subplot(2, 1, 2);
    hold on
    grid on

    plot1 = plot(T, Y(:, 2)); 
    hold on
    plot2 = plot(T, Y(:, 3)); 
    hold on
    plot3 = plot(T, Y(:, 4));
     
    title('Variation of n, m, h using HH model for I ext = 0 uA/cm2 starting near Eqbm position');
    xlabel('Time (in ms)');     
    ylabel('Ion-Channel Opening Probabilities');
    ylim([0 1]);
    legend([plot1, plot2, plot3], ["n(t)", "m(t)", "h(t)"]);
    hold off
   
    % Determining threshold Current Impulse that produces a spike

    I_impulse = 0:0.1:15;
    I_ext = 0;
    temp_arr = zeros(numel(I_impulse), 2);
    T_span = [0 100];

    for i = 1:numel(I_impulse)
    
       V_init = double(eqbm.V_mem) + (I_impulse(i)/C);
       init_condn = [V_init, double(eqbm.n), double(eqbm.m), double(eqbm.h)];
       [~, Y] =  HH_equations_for_given_current(I_ext, T_span, init_condn);

       temp_arr(i, 1) = I_impulse(i);
       temp_arr(i, 2) = max(Y(:, 1));

    end

    figure(13)
    hold on
    grid on

    plot(temp_arr(:, 1), temp_arr(:, 2));
    title('Effect of depolarizing Current Impulses from V rest = -60 mV in HH model (Threshold Behaviour)');
    xlabel('Depolarizing Current Impulse (in uA/cm2)');     
    ylabel('Max Membrane Potential (in mV)');
    hold off

    % Question 15: Model stabiity at Equilibrium for steady I_ext in range 8-12 uA/cm2 in HH model 
    fprintf("\nQUESTION 15:\n");
   
    I_ext = 8:12;
    syms n m h
    T_span = [0 300];

    figure(14)

    for i = 1:numel(I_ext)

        f1 = (1/C)*(I_ext(i) + (-G_K_const*(n^4)*(V_mem - E_K)) + (-G_Na_const*(m^3)*h*(V_mem - E_Na)) + (-G_L_const*(V_mem - E_L_rest)));         % dV/dt
        f2 = alpha_n*(1 - n) - beta_n*n;                                                                                                          % dn/dt
        f3 = alpha_m*(1 - m) - beta_m*m;                                                                                                          % dm/dt
        f4 = alpha_h*(1 - h) - beta_h*h;                                                                                                          % dh/dt
        eqbm = vpasolve([f1 == 0, f2 == 0, f3 == 0, f4 == 0], [V_mem, n, m, h]);
    
        fprintf("\nThe Equilibrium condition for HH model for steady I_ext = %f uA/cm2:\n", I_ext(i));
        fprintf("V_eq = %f mV, n_eq = %f, m_eq = %f, h_eq = %f\n", eqbm.V_mem, eqbm.n, eqbm.m, eqbm.h);
        
        % Checking the stability at equilibrium by starrting slighly away from
        % the Equilibrium position
        init_condn = [double(eqbm.V_mem) + 0.1, double(eqbm.n) + 0.01, double(eqbm.m) + 0.001, double(eqbm.h) + 0.01];
        [T, Y] = HH_equations_for_given_current(I_ext(i), T_span, init_condn);
    
        % Action Potential plot
 
        subplot(3, 2, i);
        hold on
        grid on
        plot(T, Y(:, 1));
        
        str_temp = sprintf('Stability Check for steady I external = %f uA/cm^2', I_ext(i));
        title(str_temp);
        xlabel('Time (in ms)');     
        ylabel('Membrane Potential (in mV)');
    
        hold off
    end

    % Question 16:  Action Potential for the V-n reduced HH model
    fprintf("\nQUESTION 16:\n");
    
    I_ext = 10.0;
    T_span = [0 80];

    V_init = -60;
    n_init = double(subs(n_inf_V_rest, {sym('V_mem')}, {V_init}));

    init_condn = [V_init + 0.1, n_init + 0.001];
    
    [T, Y] = V_n_reduced_HH_equations_for_given_current(I_ext, T_span, init_condn);

    figure(15)
    hold on
    grid on

    subplot(2, 1, 1);
    grid on
    hold on
    plot(T, Y(:, 1));

    title('Action Potentials using V-n reduced HH model for I ext = 10 uA/cm2 starting near Eqbm position');
    xlabel('Time (in ms)');     
    ylabel('Membrane Potential (in mV)');
    hold off

    % Plot showing variation of n, m, h with time corresponding to above
    % Action Potentials for V-n reduced model
    
    subplot(2, 1, 2);
    hold on 
    grid on
    
    % plot of n as a function of t
    plot1 = plot(T, Y(:, 2));   
    hold on
 
    % plot for m as a function of t
    alpha_m = -0.1*(V_mem + 35)/(exp(-(V_mem + 35)/10) - 1);
    beta_m = 4*(exp(-(V_mem + 60)/18));
    m = alpha_m/(alpha_m + beta_m);

    m_arr = zeros(numel(Y(:, 1)), 1);

    for i = 1:numel(Y(:, 1))
        m_arr(i) = subs(m, {sym('V_mem')}, {Y(i, 1)});
    end
   
    plot2 = plot(T, m_arr);

    % plot for h as a function of t
    h_arr = zeros(numel(Y(:, 2)), 1);

    for i = 1:numel(Y(:, 2))
        h_arr(i) = 1 - Y(i, 2);
    end
    
    hold on
    plot3 = plot(T, h_arr);

    title('Variation of n, m, h using V-n reduced HH model for I ext = 10 uA/cm2 starting near Eqbm position');
    xlabel('Time (in ms)');     
    ylabel('Ion-Channel Opening Probabilities');
    ylim([0 1]);
    legend([plot1, plot2, plot3], ["n(t)", "m(t)", "h(t)"]);
    hold off

    % Effect of depolarizing current impulse in V-n reduced HH model (Thresold Behaviour)

    I_impulse = 0:0.1:15;
    T_span = [0 100];
    I_ext = 0;

    temp_arr = zeros(numel(I_impulse), 2);

    for i = 1:numel(I_impulse)

        V_init = double(eqbm.V_mem) + (I_impulse(i)/C);
        init_condn = [V_init, double(eqbm.n)];
        [~, Y] =  V_n_reduced_HH_equations_for_given_current(I_ext, T_span, init_condn);
        
        temp_arr(i, 1) = I_impulse(i);
        temp_arr(i, 2) = max(Y(:, 1));

    end

    figure(16)
    hold on 
    grid on

    plot(temp_arr(:, 1), temp_arr(:, 2));
    title('Effect of depolarizing Current Impulses from V rest = -60 mV in V-n reduced HH model (Threshold Behaviour)');
    xlabel('Depolarizing Current Impulse (in uA/cm2)');     
    ylabel('Max Membrane Potential (in mV)');

    hold off
    
    % Model stabiity at Equilibrium for steady I_ext in range 8-12 uA/cm2 in V-n reduced HH model 

    I_ext = 8:12;
    syms n 
    T_span = [0 50];

    figure(17)

    fprintf("\n V-n REDUCED HH SYSTEM OF EQUATIONS \n\n")

    for i = 1:numel(I_ext)

        alpha_m = -0.1*(V_mem + 35)/(exp(-(V_mem + 35)/10) - 1);
        beta_m = 4*(exp(-(V_mem + 60)/18));
        m_inf_V = alpha_m/(alpha_m + beta_m);

        f1 = (1/C)*(I_ext(i) + (-G_K_const*(n^4)*(V_mem - E_K)) + (-G_Na_const*(m_inf_V^3)*(1-n)*(V_mem - E_Na)) + (-G_L_const*(V_mem - E_L_rest)));        % dV/dt
        f2 = alpha_n*(1 - n) - beta_n*n;                                                                                                                    % dn/dt
                                                                                                      
        eqbm = vpasolve([f1 == 0, f2 == 0], [V_mem, n]);
    
        fprintf("\nThe Equilibrium condition for V-n reduced system HH model for steady I_ext = %f uA/cm2:\n", I_ext(i));
        fprintf("V_eq = %f mV, n_eq = %f\n", eqbm.V_mem, eqbm.n);
        
        % Checking the stability at equilibrium by starrting slighly away from
        % the Equilibrium position
        init_condn = [double(eqbm.V_mem) + 0.1, double(eqbm.n) + 0.01];
        [T, Y] = V_n_reduced_HH_equations_for_given_current(I_ext(i), T_span, init_condn);
    
        % Action Potential plot
 
        subplot(3, 2, i);
        hold on 
        grid on

        plot(T, Y(:, 1));
        
        str_temp = sprintf('Stability Check for steady I external = %f uA/cm^2', I_ext(i));
        title(str_temp);
        xlabel('Time (in ms)');     
        ylabel('Membrane Potential (in mV)');
        hold off

    end
  
    % Question 17: Anode Break Excitation phenomemon
 
    T_span = [0 50];
    V_init = -60;
    n_init = double(subs(n_inf_V_rest, {sym('V_mem')}, {V_init}));
    m_init = double(subs(m_inf_V_rest, {sym('V_mem')}, {V_init}));
    h_init = double(subs(h_inf_V_rest, {sym('V_mem')}, {V_init}));

    init_condn = [V_init, n_init, m_init, h_init];

    [T, Y] = ode15s(@HH_model_anode_break, T_span, init_condn);

    % Action Potential plot
    figure(18)
    hold on 
    grid on

    plot(T, Y(:, 1));
    title('Action Potentials using HH model for I ext = -3 uA/cm2 for 20ms to demonstrate Anode Break Excitation');
    xlabel('Time (in ms)');     
    ylabel('Membrane Potential (in mV)');
    xline([10 30], LineStyle="-", Label= {'Onset of An Inward Current for 20 ms', 'Inward Current Injection Ends'}, LabelVerticalAlignment="middle");
   
    hold off    

    % Question 18:  Anode Break Excitation in V-m reduced system
    fprintf("\nQUESTION 18:\n");
    fprintf("\n V-m REDUCED HH SYSTEM OF EQUATIONS \n\n")

    syms V_mem

    alpha_n = -0.01*(V_mem + 50)/(exp(-(V_mem + 50)/10) - 1);
    beta_n = 0.125*exp(-(V_mem + 60)/80);
    n_inf_V_rest = alpha_n/(alpha_n + beta_n);
    
    alpha_h = 0.07*exp(-(V_mem + 60)/20);
    beta_h = 1/(exp(-(V_mem + 30)/10) + 1);
    h_inf_V_rest = alpha_h/(alpha_h + beta_h);
    
    % Case 1: V = V_rest = -60 mV, I_ext = 0 uA/cm2
    
    fprintf("\nCASE 1 :\n\n");
    % Plotting nullclines
    V_init = -60;
    I_ext = 0.0;
    n = double(subs(n_inf_V_rest, {sym('V_mem')}, {V_init}));
    h = double(subs(h_inf_V_rest, {sym('V_mem')}, {V_init}));

    v_temp = linspace(-80,60, 400);

    if v_temp == -35
        alpha_m = 1.0;
    else
         alpha_m = -0.1*(v_temp + 35)./(exp(-(v_temp + 35)/10) - 1);
    end

    beta_m = 4*(exp(-(v_temp + 60)./18));
    m_nullcline = (alpha_m)./(alpha_m + beta_m);
    V_nullcline = @(V) (((I_ext + (-G_K_const*(n^4)*(V - E_K)) + (-G_L_const*(V - E_L_rest)))./(G_Na_const*h*(V - E_Na)))^(1/3));

    figure(19)
    hold on 
    grid on

    plot1 = plot(v_temp, m_nullcline, 'Color', 'red', 'LineStyle','-', LineWidth =1.5);
    hold on
    plot2 = fplot(@(V) V_nullcline(V), [-80 60], 'Color', 'blue', 'LineStyle','-', 'LineWidth',1.5);
    xlim([-80, 60]);    ylim([0, 1.1]);

    % Finding Equilibrium points

    eqbm_pts = zeros(3, 2);    start_pt = zeros(3, 2);

    start_pt(1, :) = [-65.0, 0.025];
    start_pt(2, :) = [-52.0, 0.1];
    start_pt(3, :) = [48.0, 0.9];
                                                                                                        
    for i = 1:3

        syms V_mem m_temp

        alpha_m = -0.1*(V_mem + 35)/(exp(-(V_mem + 35)/10) - 1);
        beta_m = 4*(exp(-(V_mem + 60)/18));
        f1 = (I_ext + (-G_K_const*(n^4)*(V_mem - E_K)) + (-G_Na_const*(m_temp^3)*h*(V_mem - E_Na)) + (-G_L_const*(V_mem - E_L_rest)));                      % dV/dt
        f2 = alpha_m*(1 - m_temp) - beta_m*m_temp;                                                                                                          % dm/dt

        eqbm = vpasolve([f1 == 0, f2 == 0], [V_mem, m_temp], [start_pt(i, 1), start_pt(i, 2)]);
    
        fprintf("\nEquilibrium point %d for V-m reduced HH model (Case 1) :\n", i);
        fprintf("V_eq = %f mV    m_eq = %f\n", eqbm.V_mem, eqbm.m_temp);

        eqbm_pts(i, 1) = eqbm.V_mem;
        eqbm_pts(i, 2) = eqbm.m_temp;
        
        % Stability of Equilibrium
    
        % f1 = dv/dt; f2 = dm/dt
        f1 = (1/C)*(I_ext + (-G_K_const*(n^4)*(V_mem - E_K) + (-G_Na_const*(m_temp^3)*h*(V_mem- E_Na) )) + (-G_L_const*(V_mem - E_L_rest)));
        f2 = alpha_m*(1- m_temp) - beta_m*m_temp;
    
        j_mat = jacobian([f1; f2], [V_mem; m_temp]);
        j_mat_at_eqlbrm = subs(j_mat, {sym('V_mem'), sym('m_temp')}, {eqbm_pts(i, 1), eqbm_pts(i, 2)});
    
        eigen_vals = double(eig(j_mat_at_eqlbrm));
        res = check_stability(eigen_vals(1), eigen_vals(2));
        fprintf("The given equilibrium point is a %s.\n", res);

    end
    
    plot3 = plot(eqbm_pts(1, 1), eqbm_pts(1, 2), 'Marker','o', 'Color', 'black', 'MarkerFaceColor','green');
    plot4 = plot(eqbm_pts(2, 1), eqbm_pts(2, 2), 'Marker','o', 'Color', 'black', 'MarkerFaceColor','cyan');
    plot5 = plot(eqbm_pts(3, 1), eqbm_pts(3, 2), 'Marker','o', 'Color', 'black', 'MarkerFaceColor','magenta');
   
    % Case 2: 
    fprintf("\nCASE 2 :\n\n");

    I_ext = -3.0;
    syms V_mem n_temp h_temp m_temp

    f1 = (1/C)*(I_ext + (-G_K_const*(n_temp^4)*(V_mem - E_K)) + (-G_Na_const*(m_temp^3)*h_temp*(V_mem - E_Na)) + (-G_L_const*(V_mem - E_L_rest)));      % dV/dt
    f2 = alpha_n*(1 - n_temp) - beta_n*n_temp;                                                                                                          % dn/dt
    f3 = alpha_m*(1 - m_temp) - beta_m*m_temp;                                                                                                          % dm/dt
    f4 = alpha_h*(1 - h_temp) - beta_h*h_temp;                                                                                                          % dh/dt
    eqbm = vpasolve([f1 == 0, f2 == 0, f3 == 0, f4 == 0], [V_mem, n_temp , m_temp, h_temp]);
  
    V_init = eqbm.V_mem;
    n = double(subs(n_inf_V_rest, {sym('V_mem')}, {V_init}));
    h = double(subs(h_inf_V_rest, {sym('V_mem')}, {V_init}));

    % Plotting nullclines  
   
    V_nullcline = @(V) (((I_ext + (-G_K_const*(n^4)*(V - E_K)) + (-G_L_const*(V - E_L_rest)))./(G_Na_const*h*(V - E_Na)))^(1/3));
    
    hold on
    plot6 = fplot(@(V) V_nullcline(V), [-80 60], 'LineStyle','-', 'LineWidth',1.5);
 
    % Finding Equilibrium points

    eqbm_pts = zeros(3, 2);    start_pt = zeros(3, 2);

    start_pt(1, :) = [-65.0, 0.025];
    start_pt(2, :) = [-52.0, 0.1];
    start_pt(3, :) = [48.0, 0.9];
                                                                                                        
    for i = 1:3

        syms V_mem m_temp

        alpha_m = -0.1*(V_mem + 35)/(exp(-(V_mem + 35)/10) - 1);
        beta_m = 4*(exp(-(V_mem + 60)/18));
        f1 = (I_ext + (-G_K_const*(n^4)*(V_mem - E_K)) + (-G_Na_const*(m_temp^3)*h*(V_mem - E_Na)) + (-G_L_const*(V_mem - E_L_rest)));                      % dV/dt
        f2 = alpha_m*(1 - m_temp) - beta_m*m_temp;                                                                                                          % dm/dt

        eqbm = vpasolve([f1 == 0, f2 == 0], [V_mem, m_temp], [start_pt(i, 1), start_pt(i, 2)]);
    
        fprintf("\nEquilibrium point %d for V-m reduced HH model (Case 2) :\n", i);
        fprintf("V_eq = %f mV    m_eq = %f\n", eqbm.V_mem, eqbm.m_temp);

        eqbm_pts(i, 1) = eqbm.V_mem;
        eqbm_pts(i, 2) = eqbm.m_temp;
        
        % Stability of Equilibrium
    
        % f1 = dv/dt; f2 = dm/dt
        f1 = (1/C)*(I_ext + (-G_K_const*(n^4)*(V_mem - E_K) + (-G_Na_const*(m_temp^3)*h*(V_mem- E_Na) )) + (-G_L_const*(V_mem - E_L_rest)));
        f2 = alpha_m*(1- m_temp) - beta_m*m_temp;
  
        j_mat = jacobian([f1; f2], [V_mem; m_temp]);
        j_mat_at_eqlbrm = subs(j_mat, {sym('V_mem'), sym('m_temp')}, {eqbm_pts(i, 1), eqbm_pts(i, 2)});
    
        eigen_vals = double(eig(j_mat_at_eqlbrm));
        res = check_stability(eigen_vals(1), eigen_vals(2));
        fprintf("The given equilibrium point is a %s.\n", res);

    end
    
    hold on
    plot7 = plot(eqbm_pts(1, 1), eqbm_pts(1, 2), 'Marker','o', 'Color', 'black', 'MarkerFaceColor','blue');
    plot8 = plot(eqbm_pts(2, 1), eqbm_pts(2, 2), 'Marker','o', 'Color', 'black', 'MarkerFaceColor','red');
    plot9 = plot(eqbm_pts(3, 1), eqbm_pts(3, 2), 'Marker','o', 'Color', 'black', 'MarkerFaceColor', 'yellow');

    title('Nullclines and Equilibrium points for V-m reduced HH model (Case 1 and 2)');
    xlabel('V');    ylabel('m');

    legend([plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9], ["m-nullcline (same for C1 and C2)", "V-nullcline (C1)", "Equilibrium Point 1 - Stable (C1)", "Equilibrium Point 2 - Saddle (C1)", "Equilibrium Point 3 - Stable(C1)", "V-nullcline (C2)", "Equilibrium Point 1 - Stable (C2)", "Equilibrium Point 2 - Saddle (C2)", "Equilibrium Point 3 - Stable(C2)"], "Location","north");    
    hold off


    fprintf("\n ------------------------ END OF PROJECT -----------------------------\n")

end


function derivative_vec = MLE_phi_1(t, y)

    g_Ca = 4.4; g_K = 8.0; g_L = 2.0;
    V_Ca = 120; V_K = -84; V_L = -60; phi_1 = 0.02;
    V1 = -1.2; V2 = 18; V3 = 2; V4 = 30; C = 20; I_ext = 0.0;

    derivative_vec = zeros(2, 1);
    % First component of derivative_vec => dV/dt
    % Second component of derivative_vec => dW/dt
    derivative_vec(1) = (1./C)*(I_ext - g_Ca*(0.5*(1 + tanh((y(1) - V1)./V2)))*(y(1) - V_Ca) - g_K*y(2)*(y(1) - V_K) - g_L*(y(1) - V_L));
    derivative_vec(2) = (phi_1*(0.5*(1 + tanh((y(1) - V3)./V4)) - y(2)))./(1./cosh((y(1) - V3)./(2*V4)));

end

function derivative_vec = MLE_phi_2(t, y)

    g_Ca = 4.4; g_K = 8.0; g_L = 2.0;
    V_Ca = 120; V_K = -84; V_L = -60; phi_2 = 0.04;
    V1 = -1.2; V2 = 18; V3 = 2; V4 = 30; C = 20; I_ext = 0.0;

    derivative_vec = zeros(2, 1);
    derivative_vec(1) = (1./C)*(I_ext - g_Ca*(0.5*(1 + tanh((y(1) - V1)./V2))).*(y(1) - V_Ca) - g_K*y(2)*(y(1) - V_K) - g_L*(y(1) - V_L));
    derivative_vec(2) = (phi_2*(0.5*(1 + tanh((y(1) - V3)./V4)) - y(2)))./(1./cosh((y(1) - V3)./(2*V4)));

end

function derivative_vec = MLE_phi_3(t, y)

    g_Ca = 4.4; g_K = 8.0; g_L = 2.0;
    V_Ca = 120; V_K = -84; V_L = -60; phi_3 = 0.01;
    V1 = -1.2; V2 = 18; V3 = 2; V4 = 30; C = 20; I_ext = 0.0;

    derivative_vec = zeros(2, 1);
    derivative_vec(1) = (1./C)*(I_ext - g_Ca*(0.5*(1 + tanh((y(1) - V1)./V2)))*(y(1) - V_Ca) - g_K*y(2)*(y(1) - V_K) - g_L*(y(1) - V_L));
    derivative_vec(2) = (phi_3*(0.5*(1 + tanh((y(1) - V3)./V4)) - y(2)))./(1./cosh((y(1) - V3)./(2*V4)));

end

function derivative_vec = MLE_steady_current(t, y)

    g_Ca = 4.4; g_K = 8.0; g_L = 2.0;
    V_Ca = 120; V_K = -84; V_L = -60; phi = 0.02;
    V1 = -1.2; V2 = 18; V3 = 2; V4 = 30; C = 20; I_ext = 86.0;

    derivative_vec = zeros(2, 1);
    derivative_vec(1) = (1./C)*(I_ext - g_Ca*(0.5*(1 + tanh((y(1) - V1)./V2)))*(y(1) - V_Ca) - g_K*y(2)*(y(1) - V_K) - g_L*(y(1) - V_L));
    derivative_vec(2) = (phi*(0.5*(1 + tanh((y(1) - V3)./V4)) - y(2)))./(1./cosh((y(1) - V3)./(2*V4)));

end

function derivative_vec = MLE_steady_current_time_reverse(t, y)

    g_Ca = 4.4; g_K = 8.0; g_L = 2.0;
    V_Ca = 120; V_K = -84; V_L = -60; phi = 0.02;
    V1 = -1.2; V2 = 18; V3 = 2; V4 = 30; C = 20; I_ext = 86.0;

    derivative_vec = zeros(2, 1);
    derivative_vec(1) = (-1./C)*(I_ext - g_Ca*(0.5*(1 + tanh((y(1) - V1)./V2)))*(y(1) - V_Ca) - g_K*y(2)*(y(1) - V_K) - g_L*(y(1) - V_L));
    derivative_vec(2) = (-1*phi*(0.5*(1 + tanh((y(1) - V3)./V4)) - y(2)))./(1./cosh((y(1) - V3)./(2*V4)));

end


function [T, Y] = MLE_for_given_steady_I_param_set1(i_ext, V_start, W_start)

    T_span = [0 2000];
    
    function derivative_vec = MLE_steady_current_temp(t, y)

        g_Ca = 4.4; g_K = 8.0; g_L = 2.0;
        V_Ca = 120; V_K = -84; V_L = -60; phi = 0.02;
        V1 = -1.2; V2 = 18; V3 = 2; V4 = 30; C = 20;
    
        derivative_vec = zeros(2, 1);
        derivative_vec(1) = (1./C)*(i_ext - g_Ca*(0.5*(1 + tanh((y(1) - V1)./V2)))*(y(1) - V_Ca) - g_K*y(2)*(y(1) - V_K) - g_L*(y(1) - V_L));
        derivative_vec(2) = (phi*(0.5*(1 + tanh((y(1) - V3)./V4)) - y(2)))./(1./cosh((y(1) - V3)./(2*V4)));

    end
    
    [T, Y] = ode15s(@MLE_steady_current_temp, T_span, [V_start, W_start]);

end

function [T, Y] = MLE_for_given_steady_I_param_set2(i_ext, V_start, W_start, T_span)
  

    function derivative_vec = MLE_steady_current_temp(t, y)

        g_Ca = 4.0; g_K = 8.0; g_L = 2.0;
        V_Ca = 120; V_K = -84; V_L = -60; phi = 0.0667;
        V1 = -1.2; V2 = 18; V3 = 12; V4 = 17.4; C = 20;
    
        derivative_vec = zeros(2, 1);
        derivative_vec(1) = (1./C)*(i_ext - g_Ca*(0.5*(1 + tanh((y(1) - V1)./V2)))*(y(1) - V_Ca) - g_K*y(2)*(y(1) - V_K) - g_L*(y(1) - V_L));
        derivative_vec(2) = (phi*(0.5*(1 + tanh((y(1) - V3)./V4)) - y(2)))./(1./cosh((y(1) - V3)./(2*V4)));

    end
    
    [T, Y] = ode15s(@MLE_steady_current_temp, T_span, [V_start, W_start]);

end

function rate = calculate_firing_rate(T, Y)

    V_threshold_for_AP = 10.0;
    
    count = 0;

    for i = 1: numel(Y(:, 1))-1
        if (Y(i, 1) < V_threshold_for_AP && Y(i+1,1) >= V_threshold_for_AP)
            count = count + 1;
        elseif (Y(i, 1) > V_threshold_for_AP && Y(i+1,1) <= V_threshold_for_AP)
            count = count + 1;
        end
    end

    count = count/2;
    rate = 1000*(count/(max(T) - min(T)));   % multiply by 1000 as time in denominator in milliseconds
end

% HH MODEL FUNCTIONS

function [T, Y] = HH_equations_for_given_current(I_ext, T_span, Init_condn)
    
    function deriv_vec = HH_model(t, y)
        
        % y(1) = V; y(2) = n; y(3) = m; y(4) = h;

        G_K_const = 36;     G_Na_const = 120;    G_L_const = 0.3;
        E_K = -72;          E_Na = 55;           E_L = -49.4;
        C = 1;              epsilon = 0.000001;

        alpha_n = (-0.01*(y(1) + 50) + epsilon)/((exp(-(y(1) + 50)/10) - 1) + epsilon);
        beta_n = 0.125*exp(-(y(1) + 60)/80);
        
        alpha_m = (-0.1*(y(1) + 35) + epsilon)/((exp(-(y(1) + 35)/10) - 1) + epsilon);
        beta_m = 4*(exp(-(y(1) + 60)/18));
        
        alpha_h = 0.07*exp(-(y(1) + 60)/20);
        beta_h = 1/(exp(-(y(1) + 30)/10) + 1);

        deriv_vec = zeros(4, 1);

        deriv_vec(1) = (1/C)*(I_ext + (-G_K_const*(y(2)^4)*(y(1) - E_K)) + (-G_Na_const*(y(3)^3)*y(4)*(y(1) - E_Na)) + (-G_L_const*(y(1) - E_L)));   % dV/dt
        deriv_vec(2) = alpha_n*(1 - y(2)) - beta_n*y(2);                                                                                             % dn/dt
        deriv_vec(3) = alpha_m*(1 - y(3)) - beta_m*y(3);                                                                                             % dm/dt
        deriv_vec(4) = alpha_h*(1 - y(4)) - beta_h*y(4);                                                                                             % dh/dt
    end

    [T, Y] = ode15s(@HH_model, T_span, Init_condn);

end

function [T, Y] = V_n_reduced_HH_equations_for_given_current(I_ext, T_span, Init_condn)
    
    function deriv_vec = V_n_reduced_HH_model(t, y)
        
        % y(1) = V; y(2) = n;

        G_K_const = 36;     G_Na_const = 120;    G_L_const = 0.3;
        E_K = -72;          E_Na = 55;           E_L = -49.4;
        C = 1;              epsilon = 0.000001;

        alpha_n = (-0.01*(y(1) + 50) + epsilon)/((exp(-(y(1) + 50)/10) - 1) + epsilon);
        beta_n = 0.125*exp(-(y(1) + 60)/80);
        
        alpha_m = (-0.1*(y(1) + 35) + epsilon)/((exp(-(y(1) + 35)/10) - 1) + epsilon);
        beta_m = 4*(exp(-(y(1) + 60)/18));
        m_inf_V = (alpha_m)/(alpha_m + beta_m);
    

        deriv_vec = zeros(2, 1);

        deriv_vec(1) = (1/C)*(I_ext + (-G_K_const*(y(2)^4)*(y(1) - E_K)) + (-G_Na_const*(m_inf_V^3)*(1 - y(2))*(y(1) - E_Na)) + (-G_L_const*(y(1) - E_L)));   % dV/dt
        deriv_vec(2) = alpha_n*(1 - y(2)) - beta_n*y(2);                                                                                                    % dn/dt
    end

    [T, Y] = ode15s(@V_n_reduced_HH_model, T_span, Init_condn);

end

function deriv_vec = HH_model_anode_break(t, y)

    G_K_const = 36;     G_Na_const = 120;    G_L_const = 0.3;
    E_K = -72;          E_Na = 55;           E_L = -49.4;
    C = 1;              epsilon = 0.000001;

    alpha_n = (-0.01*(y(1) + 50) + epsilon)/((exp(-(y(1) + 50)/10) - 1) + epsilon);
    beta_n = 0.125*exp(-(y(1) + 60)/80);
    
    alpha_m = (-0.1*(y(1) + 35) + epsilon)/((exp(-(y(1) + 35)/10) - 1) + epsilon);
    beta_m = 4*(exp(-(y(1) + 60)/18));
    
    alpha_h = 0.07*exp(-(y(1) + 60)/20);
    beta_h = 1/(exp(-(y(1) + 30)/10) + 1);

    deriv_vec = zeros(4, 1);

    if (t < 10.0 || t > 30.0)
        I_ext = 0;
    else
        I_ext = -3;
    end

    deriv_vec(1) = (1/C)*(I_ext + (-G_K_const*(y(2)^4)*(y(1) - E_K)) + (-G_Na_const*(y(3)^3)*y(4)*(y(1) - E_Na)) + (-G_L_const*(y(1) - E_L)));   % dV/dt
    deriv_vec(2) = alpha_n*(1 - y(2)) - beta_n*y(2);                                                                                             % dn/dt
    deriv_vec(3) = alpha_m*(1 - y(3)) - beta_m*y(3);                                                                                             % dm/dt
    deriv_vec(4) = alpha_h*(1 - y(4)) - beta_h*y(4);                                                                                             % dh/dt 

end

% function to check stability of equilibrium point
function result = check_stability(eigenval1, eigenval2)

    result = "";
    
    if (imag(eigenval1) == 0)
        
        if (eigenval1 < 0.0 && eigenval2 < 0.0)
            result = "Stable Equilibrium";   
        elseif (eigenval1 > 0.0 && eigenval2 > 0.0)
            result = "Unstable Equilibrium";
        elseif (eigenval1*eigenval2 < 0.0)
            result = "Saddle Point";
        end

    else

        if (real(eigenval1) < 0.0)
            result = "Stable Spiral";
        else
            result = "Unstable Spiral";
        end

    end
end
