% Pranav Nyati
% 20CS30037
% Proj-3


function PRANAV_20CS30037_PROJECT3()


    Vars = {'All_Spike_Times', 'Stimulus'};
    
    data_mat = load("data_cn_project_iii_a22.mat", Vars{:});
    
    spike_times = data_mat.All_Spike_Times;
    stimulus = data_mat.Stimulus;
    
    exptn_stimuli = mean(stimulus);
    fprintf("The expectation of the stimulus over time : %.6f\n", exptn_stimuli);
    
    var_stimuli = var(stimulus);
    fprintf("The variance of the stimulus over time : %.6f\n", var_stimuli);
    
    num_repeat = 50;
    
    % Question 1: Stimulus nature

    time_bins = 0:0.001:19.999;
    figure(1)
    plot(time_bins, stimulus);
    title('Stimulus as a function of time for 0-20s');
    xlabel('Time in seconds');
    ylabel('Stimulus (in dB)');
    grid on;
    hold off
    
    % Autocorrelation function:
    tau_vals = -0.050: 0.001: 0.050;
    auto_corr_tau = zeros(101, 1);
    
    
    for i = 1:numel(tau_vals)
        for j = 1:numel(time_bins)
            if(j + 1000*tau_vals(i) > numel(time_bins) || j + 1000*tau_vals(i) < 1)
                continue;
            else
                auto_corr_tau(i) = auto_corr_tau(i) + (0.001*stimulus(j)*stimulus(round(j + 1000*tau_vals(i))));
            end
        end
        
        % dividing by the total time period of the stimulus
        auto_corr_tau(i) = auto_corr_tau(i)/20;
    end
    
    figure(2)
    plot(tau_vals*1000, auto_corr_tau);
    title('Plot of Autocorrelation of the Stimulus with varying Tau');
    xlabel('Tau (in milliseconds)');
    ylabel('Autocorrelation value');
    grid on;
    hold off
    
    
    % Question 2: PSTH and mean firing rate
    
    time_bins = 0:0.001:20;
    spike_counts = zeros(4, numel(time_bins) - 1);
    
    tot_spike_count = 0;
    
    for i = 1:4
        for j = 1:num_repeat
            temp = spike_times{i, j};
            for k = 1:numel(temp)
                tot_spike_count = tot_spike_count + numel(temp);
                index = ceil(1000*temp(1, k));
                spike_counts(i, index) = spike_counts(i, index) + 1; 
            end
            
        end
    
    end
    
    spike_rates = zeros(4, numel(time_bins) - 1);
    
    for i = 1:size(spike_rates, 1)
        for j = 1:size(spike_rates, 2)
            spike_rates(i, j) = (spike_counts(i, j)/num_repeat)*1000;
        end
    end
    

    % Plots of PSTH for the 4 neurons
    figure(3)
   
    time_steps = 0.001:0.001:20;
    
    PSTH_neuron_1 = plot(time_steps, spike_rates(1, :));
    ylim([0 160]);
    title('PSTH: neuron 1');
    xlabel('Time (in seconds)');
    ylabel('Rate (spikes/second)');
  
    grid on;
    hold off
    
    figure(4)

    time_steps = 0.001:0.001:20;
    
    PSTH_neuron_2 = plot(time_steps, spike_rates(2, :));
    ylim([0 160]);    
    title('PSTH: neuron 2');
    xlabel('Time (in seconds)');
    ylabel('Rate (spikes/second)');
    grid on;
    hold off
    
    figure(5)
   
    time_steps = 0.001:0.001:20;
    
    PSTH_neuron_3 = plot(time_steps, spike_rates(3, :));
    title('PSTH: neuron 3');
    xlabel('Time (in seconds)');
    ylabel('Rate (spikes/second)');
    ylim([0 160]);
    grid on;
    hold off
    
    figure(6)
    grid on;
    
    time_steps = 0.001:0.001:20;
    
    PSTH_neuron_4 = plot(time_steps, spike_rates(4, :));
    title('PSTH: neuron 4');
    xlabel('Time (in seconds)');
    ylabel('Rate (spikes/second)');
    ylim([0 160]);
    grid on;
    hold off
    
    % Question 3: Poisson or non-Poisson
    
    time_bin_size = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5];
    num_repeat = 50;
    
    mean_mat = zeros(numel(time_bin_size), 4, num_repeat);
    variance_mat = zeros(numel(time_bin_size), 4, num_repeat);
    
    
     for i = 1:numel(time_bin_size)
    
        time_bins = 0:time_bin_size(i):20;
    
        for j = 1:4
            for k = 1:num_repeat
    
                spike_counts = zeros(numel(time_bins) - 1, 1);
                temp = spike_times{j, k};
    
                for l = 1:numel(temp)
                    index = ceil((1/time_bin_size(i))*temp(1, l));
                    spike_counts(index) = spike_counts(index) + 1; 
                end 
    
               
                mean_mat(i, j, k) = mean(spike_counts);
                variance_mat(i, j, k) = var(spike_counts);
         
            end   
        end
     end
    
    
     % Plots of mean vs variance for the spike data of the 4 neurons for
     % different bin sizes to characterize their Poisson or Non-Poisson
     % nature

     for i = 1:4
         for j = 1:numel(time_bin_size)
            
             figure(6 + 6*(i-1) + j);
             scatter(squeeze(mean_mat(j, i, :)), squeeze(variance_mat(j, i,:)));
             str = sprintf('Mean vs Variance plot for neuron %d for bin size %d milliseconds', i, 1000*time_bin_size(j));
             title(str);
             xlabel('Mean');
             ylabel('Variance');
             hold off
         end
     end
    
    
    
     % Question 4: Spike Triggered Average (And Correction for Non-Gaussianity)
    
    sta_neurons = zeros(4, 100);
    num_spikes_each  = zeros(4, 1);
    stimulus_15_sec = stimulus(1: 15000);
    
    for i = 1:4
        
        for j = 1:num_repeat
            
            % all spike times in time interval 0-15 seconds of ith neuron in jth repeation 
            temp = spike_times{i, j}(spike_times{i, j} < 15);
            
            % updating the no of spikes for a neuron over the 15 seconds
            % duration of all repeatations
            num_spikes_each(i) = num_spikes_each(i) + numel(temp);
    
            for k = 1:numel(temp)
                for l = 1:100
                    sta_neurons(i, l) = sta_neurons(i, l) + stimulus_15_sec(max(1, ceil(1000*temp(k)) - l));
    
                end
            end
        end
    
        for t = 1:100
            sta_neurons(i, t) = sta_neurons(i, t)/num_spikes_each(i);
        end
    
    end
    
    % Plots of STA for the 4 neurons
    
    figure(31)
    
    subplot(2,2,1)
    plot(sta_neurons(1,:));
    ylim([-0.25 0.25]);
    title('Spike Triggered Average for 1st neuron for 100 ms stimulus window');
    xlabel('Time (in ms) before a spike ');
    ylabel('Average stimulus value over all repeats ');
    grid on;
    
    subplot(2,2,2)
    plot(sta_neurons(2,:));
    ylim([-0.25 0.25]);
    title('Spike Triggered Average for 2nd neuron for 100 ms stimulus window')
    xlabel('Time (in ms) before a spike ');
    ylabel('Average stimulus value over all repeats ');
    grid on;
    
    subplot(2,2,3)
    plot(sta_neurons(3,:));
    ylim([-0.25 0.25]);
    title('Spike Triggered Average for 3rd neuron for 100 ms stimulus window')
    xlabel('Time (in ms) before a spike ');
    ylabel('Average stimulus value over all repeats ');
    grid on;
    
    subplot(2,2,4)
    plot(sta_neurons(4,:));
    ylim([-0.25 0.25]);
    title('Spike Triggered Average for 4th neuron for 100 ms stimulus window')
    xlabel('Time (in ms) before a spike ');
    ylabel('Average stimulus value over all repeats ');
    grid on;
    hold off
    
    
    
    % Correction term for the filter h(t) (which is prop to STA) which is
    % computed above
    
    
    auto_corr_vec = autocorr(stimulus, 99);
    
    CSS = toeplitz(auto_corr_vec);
    size(CSS)
    
    corrected_filter = CSS\transpose(sta_neurons);
    size(corrected_filter)
    
    for i = 1:4
        figure(31 + i)
        plot(linspace(1, 100), corrected_filter(:, i));
        xlabel('Time (in milliseconds)');
        ylabel('h(t)');
        str_temp = sprintf('Corrected filter for neuron %d as a function of time', i);
        title(str_temp);
        grid on;
        hold off
    end
    
    


    % Question 5: Determining the output non-linearity
    
    
    % output of the linear system = x(t)*y(t)
    for i = 1:4
        y_t(i,:) = conv(stimulus(1, 1:15000), corrected_filter(:, i));
    end
    
    
    % Mean plot of lambda(t) for binned y(t) values
    
    x_avg = zeros(4,300);
    y_avg = zeros(4,300);
    
    bin_size = 40;
    bin_idx = 1:bin_size:15000;
    for i = 1:4
        for j = 1:ceil(15000/bin_size) 
            x_avg(i,j) = mean(y_t(i,bin_idx(j):bin_idx(j)+(bin_size-1)));
            y_avg(i,j) = mean(spike_rates(i,bin_idx(j):bin_idx(j)+(bin_size-1)));
        end
    
        figure(36);
    
        subplot(2,2,i);
        scatter(x_avg(i,:),y_avg(i,:));
        xlabel('y(t)');
        ylabel('λ(t)');
        title(['Estimated y(t) vs Measured λ(t) for Neuron ' num2str(i)]);
        grid on
    end
    
    % finding the nonlinear output function from y_mean and x_mean
    
    fit_const = cell(4,1);
    fit_error_values = struct( 'sse', cell( 4, 1 ), 'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
    
    figure(37)

    for i = 1:4
    
        % Fit for Neuron i
        [x_data,y_data] = prepareCurveData(x_avg(i,:),y_avg(i,:));
        
        % Using the sigmoid curve to fit and setting the fittype options
        ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' ); opts.Algorithm = 'Levenberg-Marquardt'; opts.Display = 'Off'; opts.Robust = 'Bisquare';
        opts.StartPoint = [0.75 0.25 0.5];
        
        % Fit model to data.
        [fit_const{i}, fit_error_values(i)] = fit( x_data, y_data, ft, opts );
        
        % Plot fit with data.
        subplot(2,2,i);
        p = plot( fit_const{i}, x_data, y_data);
        set(p, 'LineWidth', 2);
        str1 = sprintf('Non-linear fit for Neuron %d', i);
        title(str1);
        xlabel('y(t)')
        ylabel('λ(t)')
        grid on
    
    end
    
    
    % Question 6: Prediction performance and pruning of filter parameters
     
    
    x_test_preds = zeros(4,5000);
    y_test = zeros(4,5000);
    
    
    % Evaluating lambda(t) from the h(t) and non-linearity estimated before
    for i = 1:4
        
        % Output of linear filter (y(t))
        pred = conv(stimulus(15001:20000), corrected_filter(:, i));
        x_test_preds(i,:) = pred(1:5000);
    
        % Output of passing y(t) through the non-linearity estimated in part(5)
        for j = 1:5000
            x_test_preds(i, j) = (fit_const{i,1}.a)/(1+exp(-fit_const{i,1}.b*(x_test_preds(i, j)-fit_const{i,1}.c)));
        end
        
        % Firing rates of the neuron estimated from binning (by plotting the
        % PSTH)
        y_test(i,:) = spike_rates(i,15001:20000);
    end
    
    
    % Plotting the Estimated and Actual firing rate
    x_mean_test_preds = zeros(4,100);
    y_mean_test = zeros(4,100);
    
    bin_size = 50;
    bin_idx = 1:bin_size:5000;
    
    % averaging the rate over 100 bins of bin-size 50 ms
    for i = 1:4
    
        for j = 1:ceil(5000/bin_size) 
            x_mean_test_preds(i,j) = mean(x_test_preds(i,bin_idx(j):bin_idx(j)+(bin_size-1)));
            y_mean_test(i,j) = mean(y_test(i,bin_idx(j):bin_idx(j)+(bin_size-1)));
        end
    
        figure(38);
        subplot(2,2,i);
        plot(bin_idx,x_mean_test_preds(i,:));
        hold on;
        plot(bin_idx,y_mean_test(i,:));
    
        str1 = sprintf('Predicted λ(t) vs Actual λ(t) for neuron %d for last 5 seconds before pruning', i);
        title(str1);
        ylabel('λ(t)');
        xlabel('Time (in milliseconds)');
    
        legend('Actual Rate','Predicted Rate');
        hold off;
        grid on;
    end
    
    
    %  Evaluating the r_squared value between actual rates and rates predicted
    %  by the model
    
    r_sq_vals = zeros(4,1);
    fprintf('\n BEFORE PRUNING:  \n');
    
    for i = 1:4
        r_value = corrcoef(y_test(i,:),x_test_preds(i,:));
        r_sq_vals(i) = r_value(2)^2;
    
        fprintf('R squared value for neuron %d is',i);
        disp(r_sq_vals(i));
    end
    
    
    % Pruning of filter parameters
    
    ctr = zeros(100,4);
    cor_coef = zeros(100,4);
    
    x_test_pruned = zeros(4,5000);
    
    % No need to prune neuron 1 and 4, as their r-squared value is very low
    % They do not lead to correct prediction of rate 
    
    % Pruning only neuron 2 and neuron 3's filter params
    
    for i = 2:3
    
    
        r_sq_current = r_sq_vals(i);
        r_sq_new = r_sq_current;
        index = 0;
        count = 0;
    
        while (r_sq_current - r_sq_new) < 0.011  % the threshold here affects the level of pruning and has been set based on experimenting with different values
    
            r_sq_current = r_sq_new;
            min = max(corrected_filter(:, i));
    
            for j = 1:100
                if abs(corrected_filter(j, i)) < min && abs(corrected_filter(j, i)) ~= 0
                    min = abs(corrected_filter(j, i));
                    index = j;
                end
            end
    
            if index == 0
                break;
            end
    
            corrected_filter(index, i) = 0;
            index = 0;
    
            pred_pruned = conv(stimulus(15001:20000),corrected_filter(:, i));
            x_test_pruned(i,:) = pred_pruned(1:5000);
    
            for j = 1:5000
                x_test_pruned(i,j) = (fit_const{i,1}.a)/(1+exp(-fit_const{i,1}.b*(x_test_pruned(i,j)-fit_const{i,1}.c)));
            end
    
            r_value = corrcoef(y_test(i,:),x_test_pruned(i,:));
            r_sq_new = r_value(2)^2;
            count = count + 1;
            ctr(count,i) = count;
            cor_coef(count,i) = r_sq_new;
        end
    
    
        figure(38+i-1);
        scatter(ctr(:,i),cor_coef(:,i),'filled');
        str1 = sprintf('Plot of prediction performance vs no of iterations for neuron %d', i);
        title(str1);
        ylabel('Prediction performance');
        xlabel('Number of iterations');
    
    end
    
    
    % Ploting the filter for each neuron after pruning
    for i = 1:4  
        figure(41);
        subplot(2,2,i);
        plot(0:99, corrected_filter(:, i), 0:99, zeros(1,100));
        xlabel('Time (in milliseconds)');
        ylabel('h(t)');
        ylim([-1 1]);
        title(['Filter h(t) after pruning of neuron ' num2str(i)]);
    end
    
    
    % Fourier Transform (FT) of the filters
    FT_filter = zeros(4,100);
    
    axis = (1:100);
    
    for i = 1:4
    
        y_shift = circshift(corrected_filter(:, i),[0,length(corrected_filter(:, i))/2]);
        FT_filter(i,:) = fft(y_shift,length(corrected_filter(:, i)))/length(corrected_filter(:, i));
        Y_fft = circshift(FT_filter(i,:),[0,length(corrected_filter(:, i))/2]);
    
        figure(42);
        subplot(2,2,i);
        plot(axis,abs(Y_fft), axis,zeros(1,100));
    
        title(['Fourier transform h(f) of the filter h(t) of Neuron ' num2str(i)]);
        xlabel('Freqeuncy ');
        ylabel('h(f)');
        
    end
    
    % Prediction of rate by the filter post pruning 
    x_test = zeros(4,5000);
    
    for i = 1:4
        
        % Output of the LTI system (y(t) = x(t)*h(t) )
        pred = conv(stimulus(15001:20000), corrected_filter(:, i));
        x_test(i,:) = pred(1:5000);
    
        % Applying the non-linearity on y(t) to get rate(t) predicted by the
        % model
        for j = 1:5000
            x_test(i,j) = (fit_const{i,1}.a)/(1+exp(-fit_const{i,1}.b*(x_test(i,j)-fit_const{i,1}.c)));
        end
    end
    
    
    % Evaluation of r-squared values again after pruning of filter parameters
    fprintf('\n AFTER PRUNING:  \n');
    
    for i = 1:4
        r_value = corrcoef(y_test(i,:),x_test(i,:));
        r_sq_vals(i) = r_value(2)^2;
    
        fprintf('R squared value for neuron %d is',i);
        disp(r_sq_vals(i));
    
    end
    
    % Plot of predicted and actual rates along the two axes
    x_mean_test = zeros(4,100);
    y_mean_test = zeros(4,100);
    
    bin_size = 50;
    bin_idx = 1:bin_size:5000;
    
    for i = 1:4
    
        for j = 1:ceil(5000/bin_size) 
            x_mean_test(i,j) = mean(x_test(i,bin_idx(j):bin_idx(j)+(bin_size-1)));
            y_mean_test(i,j) = mean(y_test(i,bin_idx(j):bin_idx(j)+(bin_size-1)));
        end
    
        figure(43);
        subplot(2,2,i);
        plot(bin_idx,y_mean_test(i,:));
        hold on;
        plot(bin_idx,x_mean_test(i,:));
    
        str1 = sprintf('Predicted λ(t) vs Actual λ(t) for neuron %d for last 5 seconds before pruning', i);
        title(str1);
        ylabel('λ(t)');
        xlabel('Time (in milliseconds)');
    
        legend('Actual Rate','Predicted Rate');
        hold off;
        grid on;
    
    end
    
    
    
    % Part 2
    
    % Victor Purpura Distance 
    
    fprintf("\n\n PART-2 : VICTOR-PURPURA SPIKE DISTANCE MEASURE \n\n");

    % weights for the time difference
    q_array = [0 0.001 0.01 0.1 1 10 100];
    
    num_repeats = 100;
    
    MI_mat = zeros(4, num_repeats, length(q_array));

    % calculation of V-P SDM
    for i = 1: num_repeats

        fprintf("Iteration No. %d  in progress ............\n", i);
    
        rand_array = randperm(19901,8);
    
        for n = 1:4
    
            spike_segments = cell(8,50);
    
            for rep_idx = 1:50
    
                spike_temp = spike_times{n,rep_idx};
    
                for j = 1:8
                    spike_segments{j,rep_idx} = spike_temp(spike_temp >= rand_array(j)/1000 & spike_temp < (rand_array(j)+100)/1000);
                end
            end
            
            for idx = 1:length(q_array)
    
                confusion_mat = zeros(8,8);
    
                for j = 1:8
    
                    for rep = 1:50
                        
                        k = calc_closest_stimuli_idx(j, rep, spike_segments, q_array(idx));
                        confusion_mat(j,k) = confusion_mat(j,k)+1;
                    end
    
                end
    
                confusion_mat = confusion_mat/50;
    
                MI_mat(n,i,idx) = MI_mat(n,i,idx) + calc_MI(confusion_mat);
            end
        end
        

        fprintf("Iteration No. %d  completed! \n\n", i);
    end
    
    % calculation of confidence interval limits of 90 % from the MI matrix
    ci_fraction = 0.9;
    MI_std = std(abs(MI_mat),0,2);
    temp = 1 - ci_fraction;
    conf_intv_limits = tinv(1-temp/2, 99)*MI_std/sqrt(99);

    temp_mat1 = MI_mat(:,1,:);
    temp_mat2 = conf_intv_limits(:,1,:);

    figure(44)

    for i = 1:4

        subplot(2,2,i);

        plot(log10(q_array), temp_mat1(i,:));
        hold on
        plot(log10(q_array), temp_mat1(i,:)-temp_mat2(i,:), 'r--');
        plot(log10(q_array), temp_mat1(i,:)+temp_mat2(i,:), 'r--');

        [~,p] = max(temp_mat1(i,:));

        plt1 = plot(log10(q_array(p)), temp_mat1(i,1,p), '+');
        set(plt1, 'linewidth', 2)
        grid on
        hold off
        str1 = sprintf('Mean MI with 90 percent confidence interval for neuron %d', i);
        title(str1);
        
        xlabel('q (in logarithmic scale)')
        ylabel('Mean MI(q) (with 90% confidence intervals)')
    end

end


% ___________________________________________________________________________________________________________________________________________________
% HELPER FUNCTIONS

% function to calculate the closest stimuli index
function closest_idx = calc_closest_stimuli_idx(seg_num,rep_idx,spike_segments,q_arr_idx)

    mean_dist = zeros(1,8);

    for i = 1:8

        for rep = 1:50

            if (rep == rep_idx && i == seg_num)
                continue
            end

            mean_dist(i) = mean_dist(i) + calc_v_p_sdm(spike_segments{seg_num,rep_idx},spike_segments{i,rep},q_arr_idx);
        end
    end
   
    [~, k] = min(mean_dist);

    closest_idx = k;

end

% function to calculate victor purpura distance between two spike segments
function res = calc_v_p_sdm( spike_seg1, spike_seg2, q_arr_idx)
    
    num_spikes_1=length(spike_seg1);
    num_spikes_2=length(spike_seg2);
    
    % Dynamic Programming solution for calculating V_P spike distance
    % measure
    if q_arr_idx == 0
       res = abs(num_spikes_1-num_spikes_2);
       return

    elseif q_arr_idx == Inf
       res = num_spikes_1 + num_spikes_2;
       return
    end
    
    scr = zeros(num_spikes_1 + 1,num_spikes_2 + 1);
    scr(:,1)=(0 : num_spikes_1)';
    scr(1,:)=(0 : num_spikes_2);

    if(num_spikes_1 && num_spikes_2)

       for i = 2:num_spikes_1 + 1
          for j= 2:num_spikes_2 + 1
              scr(i,j) = min([scr(i-1,j) + 1, scr(i,j-1)+1, scr(i-1,j-1) + q_arr_idx*abs(spike_seg1(i-1)-spike_seg2(j-1))]);
          end

       end
    end

    res=scr(num_spikes_1 + 1, num_spikes_2 + 1);

end


% function to calculate MI from confusion matrix
function MI_res = calc_MI(confusion_mat)

    MI_res=0;
    % confusion matrix has entries of the form p(y/x)
    
    % calcualtion of MI from confusion matrix
    for i = 1:size(confusion_mat,1)
        for j = 1:size(confusion_mat,2)
            if(confusion_mat(i,j) ~= 0)
                MI_res = MI_res + confusion_mat(i,j)/size(confusion_mat,1)*log2(confusion_mat(i,j)/sum(confusion_mat(:,j)/size(confusion_mat, 1)));         
            end
        end
    end
end