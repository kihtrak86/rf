%--------------------------------------------------------------------------
%       Radar Simulation - VISUAL MULTI-TARGET TEST
%--------------------------------------------------------------------------
% Authors: Miguel Carralero Lanchares and Francisco Orcha Kovacs

clear; close all; clc;

%% --- Initial Parameters and Constants ---
Pdg = 1e5;          % Peak Power (W)
Gant_dB = 30;       % Antenna Gain (dB)
f = 9e9;            % Frequency (Hz)
tau = 1e-6;         % Pulse Width (s)
Ls_dB = 5;          % System Losses (dB)
Fn_dB = 4;          % Receiver Noise Figure (dB)
Pfa = 1e-6;         % Probability of False Alarm
Pd = 0.9;           % Probability of Detection
RCS_mean_calc = 1;  % Mean RCS (m^2)
theta_bw_deg = 3.0; % Antenna Beamwidth (degrees)

% Physical constants and derived values
c = physconst('LightSpeed'); 
k = physconst('Boltzmann'); 
T0 = 290;
Gant = 10^(Gant_dB/10); 
lambda = c/f; 
Ls = 10^(Ls_dB/10);
Fn = 10^(Fn_dB/10); 
Bn = 1/tau; 
theta_bw_rad = deg2rad(theta_bw_deg);

%% • First, calculate the maximum range (Rmax) for the Pdg,
%    Gant, system losses, pulse width, SNR based on Pfa and Pd,
%    and the mean radar cross-section of the targets (RCS_mean_calc).
% • Incorporating fluctuating targets (Note: Rmax is calculated for
%   non-fluctuating RCS_mean_calc, but the actual Pd will vary due to
%   the RCS fluctuation implemented later. For a rigorous Rmax calculation
%   with fluctuating targets, SNR_minimum_required should be adjusted
%   according to the Swerling model).
if exist('shnidman', 'file') && ~isempty(which('shnidman'))
    SNR_req_dB = shnidman(Pd, Pfa); 
else
    warning('shnidman not found. Using Albersheim approximation N=1 (non-fluctuating).');
    A_alb = log(0.62 / Pfa); 
    B_alb = log(Pd / (1 - Pd));
    SNR_req_dB = A_alb + 0.12 * A_alb * B_alb + 1.7 * B_alb;
end
SNR_minimum_required = 10^(SNR_req_dB / 10);

Ntot_mean_value = k * T0 * Fn * Bn;

Rmax_num = Pdg * Gant^2 * lambda^2 * RCS_mean_calc;
Rmax_den = (4*pi)^3 * Ls * Ntot_mean_value * SNR_minimum_required;
Rmax = (Rmax_num / Rmax_den)^(1/4);
fprintf('Maximum Range (Rmax, based on non-fluctuating RCS_mean_calc): %.2f km\n', Rmax/1000);

% Minimum resolution and minimum distance based on pulse width
R_res = c * tau / 2; R_min = R_res;
fprintf('Minimum Resolution: %.2f m\nMinimum Distance: %.2f m\n', R_res, R_min);

%% • From this:
% Adjust PRF to correspond to Rmax_na (unambiguous range)
Rmax_na = (4/3) * Rmax;
if Rmax >= Rmax_na
    warning('Rmax >= Rmax_na.');
    Rmax_na = Rmax * 1.1;
end
PRI = 2 * Rmax_na / c;
PRF = 1 / PRI;
fprintf('Rmax_na: %.2f km, PRF: %.2f Hz, PRI: %.4f ms\n', Rmax_na/1000, PRF, PRI*1000);

% Knowing the minimum SNR... Smin_receiver... Prx_min... Vrx_min...
Smin_receiver = SNR_minimum_required * Ntot_mean_value;
Prx_min_antenna = Smin_receiver * Ls;
fprintf('Smin_receiver: %.2e W, Prx_min_antenna: %.2e W\n', Smin_receiver, Prx_min_antenna);

% Derive the detection threshold voltage from Pfa
Vthreshold_absolute_envelope = sqrt(Ntot_mean_value * 2 * log(1/Pfa)); 
fprintf('Absolute Envelope Voltage Threshold (Vthreshold): %.2e V\n', Vthreshold_absolute_envelope);


%% --- Dynamic Simulation Setup ---
% • For N_TARGETS mobile targets
% • Without changing Pdg, Gant, system losses, pulse width, Pfa and Pd,
%   place N_TARGETS targets at random positions (which may be outside Rmax),
%   with random RCS values centered around the mean used for Rmax calculation.

sim_time = 180; 
scan_rate_deg = 72; 
scan_rate_rad = deg2rad(scan_rate_deg);
dt = PRI;
rng('shuffle');

% --- Multi-target definitions ---
N_TARGETS = 5;
target_pos = zeros(N_TARGETS, 2); 
velocity = zeros(N_TARGETS, 2);
target_colors = lines(N_TARGETS); 

max_initial_range = Rmax_na * 1.3; 
min_initial_range_factor = 0.3; 

% Initialize each target
for k_target = 1:N_TARGETS
    start_angle_rad_k = rand * 2*pi;
    start_range_k = (R_min + Rmax_na * min_initial_range_factor) + rand * (max_initial_range - (R_min + Rmax_na * min_initial_range_factor));
    
    % Increase base speeds and randomness range
    if start_range_k > Rmax_na
        speed_k = 600 + rand*300; % Range: 600 to 900 m/s
    else
        speed_k = 400 + rand*250; % Range: 400 to 650 m/s
    end
    target_pos(k_target, :) = [start_range_k*cos(start_angle_rad_k), start_range_k*sin(start_angle_rad_k)];
    
    % Initial direction always toward radar (origin)
    if norm(target_pos(k_target, :)) > 1e-3 % Avoid division by zero if starting exactly at origin
        velocity_direction_k = -target_pos(k_target, :) / norm(target_pos(k_target, :)); 
    else % If target starts at origin, assign a random direction
        random_direction_angle_k = rand * 2*pi;
        velocity_direction_k = [cos(random_direction_angle_k), sin(random_direction_angle_k)];
    end
    velocity(k_target, :) = speed_k * velocity_direction_k;
    fprintf('Target %d RANDOM Start: R=%.1f km, Ang=%.1f deg, Vel=%.1f m/s, Dir Toward Radar\n', k_target, start_range_k/1000, rad2deg(start_angle_rad_k), speed_k);
end

time_vec = 0:dt:sim_time; 
n_steps_alloc = length(time_vec); 
target_pos_history = cell(N_TARGETS, 1); 
for k_target = 1:N_TARGETS
    target_pos_history{k_target} = NaN(n_steps_alloc, 2); 
end
detections = []; 

% Recommended: at least two display panels
main_fig = figure('Position',[50 50 1700 750]); 
display_range_limit_initial = Rmax_na*1.2; 
circ_coords_theta=linspace(0,2*pi,100);

% --- Subplot position definitions [left, bottom, width, height] ---
pos_real_plot   = [0.05, 0.40, 0.43, 0.55]; 
pos_ppi_plot    = [0.53, 0.40, 0.43, 0.55]; 
pos_ascope_plot = [0.25, 0.05, 0.50, 0.28]; 

% --- Real Position Display ---
h_ax_real=axes('Parent', main_fig, 'Position', pos_real_plot); 
title(h_ax_real,'Real Position'); 
xlabel(h_ax_real,'X (m)'); 
ylabel(h_ax_real,'Y (m)');
axis(h_ax_real,'equal'); 
grid(h_ax_real, 'on'); 
hold(h_ax_real, 'on'); 
plot(h_ax_real,0,0,'k^','MarkerSize',8,'DisplayName','Radar');
plot(h_ax_real,Rmax*cos(circ_coords_theta),Rmax*sin(circ_coords_theta),'--b','LineWidth',0.5,'DisplayName','Rmax (ref)');
plot(h_ax_real,Rmax_na*cos(circ_coords_theta),Rmax_na*sin(circ_coords_theta),'--r','LineWidth',0.5,'DisplayName','Rmax_{na}');
fill(h_ax_real,[Rmax*cos(circ_coords_theta),fliplr(Rmax_na*cos(circ_coords_theta))],[Rmax*sin(circ_coords_theta),fliplr(Rmax_na*sin(circ_coords_theta))],'y','FaceAlpha',0.1,'EdgeColor','none','DisplayName','Blind Zone');
% --- Handles for multiple trajectories and positions ---
h_target_real_lines = gobjects(N_TARGETS, 1);
h_target_curr_pos = gobjects(N_TARGETS, 1);
for k_target = 1:N_TARGETS
    h_target_real_lines(k_target)=plot(h_ax_real,NaN,NaN,'-','Color', target_colors(k_target,:),'LineWidth',1.5,'DisplayName',sprintf('Target %d', k_target));
    h_target_curr_pos(k_target)=plot(h_ax_real,NaN,NaN,'o','Color', target_colors(k_target,:),'MarkerSize',6,'MarkerFaceColor',target_colors(k_target,:));
end
legend(h_ax_real,'Location','northwest');
axis(h_ax_real,display_range_limit_initial*[-1 1 -1 1]);

% --- PPI Display ---
ax_ppi=polaraxes('Parent',main_fig, 'Position', pos_ppi_plot); 
title(ax_ppi,'Radar PPI'); 
ax_ppi.ThetaZeroLocation='top'; 
ax_ppi.ThetaDir='clockwise';
ax_ppi.RLim=[0 display_range_limit_initial];
ax_ppi.Color=[0 0.1 0.15]; 
ax_ppi.GridColor=[0.5 1 0.5]; 
ax_ppi.RColor=[0.5 1 0.5]; 
ax_ppi.ThetaColor=[0.5 1 0.5]; 
ax_ppi.RAxis.Label.Color = 'k';     
ax_ppi.ThetaAxis.Label.Color = 'k'; 
ax_ppi.Title.Color = 'k';          
hold(ax_ppi, 'on');
polarplot(ax_ppi,circ_coords_theta,Rmax*ones(size(circ_coords_theta)),'--b','LineWidth',0.5,'HandleVisibility','off');
polarplot(ax_ppi,circ_coords_theta,Rmax_na*ones(size(circ_coords_theta)),'--r','LineWidth',0.5,'HandleVisibility','off');
h_scan_line_ppi=polarplot(ax_ppi,[NaN NaN],[0 display_range_limit_initial],'Color',[0.8 1 0.8],'LineWidth',1.5,'DisplayName','Scan');
h_detections_ppi_plot=polarplot(ax_ppi,NaN,NaN,'o','Color',[1 0.8 0.8],'MarkerSize',5,'MarkerFaceColor',[1 0.5 0.5],'DisplayName','Detections');
ppi_r_axis=ax_ppi.RAxis;
ppi_r_axis.Label.String='Range (km)';
ppi_r_axis.Label.Color = 'w';
ppi_r_axis.TickLabelFormat='%g';
ppi_r_axis.TickValues=linspace(0,display_range_limit_initial,5); 
ppi_r_axis.TickLabels=string(round(ppi_r_axis.TickValues/1000));
legend(ax_ppi,'TextColor','w','Location','southoutside','Orientation','horizontal');


% --- Display Type 1: A-SCOPE ---
h_ax_ascope = axes('Parent', main_fig, 'Position', pos_ascope_plot); 
title(h_ax_ascope,'A-Scope Echo (Amplitude vs Distance)');
xlabel(h_ax_ascope,'Distance (km)'); 
ylabel(h_ax_ascope,'Envelope Amplitude (V)');
hold(h_ax_ascope, 'on');
grid(h_ax_ascope, 'on');
set(h_ax_ascope, 'Color', 'k'); 
set(h_ax_ascope, 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.4 0.4 0.4]); 
h_ax_ascope.XLabel.Color = 'k'; 
h_ax_ascope.YLabel.Color = 'k'; 
h_ax_ascope.Title.Color = 'k';  

h_ascope_plot = plot(h_ax_ascope, NaN, NaN, 'g-', 'LineWidth', 1.5); 
h_ascope_threshold = yline(h_ax_ascope, Vthreshold_absolute_envelope, '--r', 'LineWidth', 1);
set(h_ascope_threshold, 'DisplayName', 'V_{threshold}');
ascope_range_display_km = Rmax_na / 1000;
xlim(h_ax_ascope, [0 ascope_range_display_km]);
num_ascope_points = 500;
ascope_range_bins_km = linspace(0, ascope_range_display_km, num_ascope_points);
ascope_base_noise_amplitude = sqrt(Ntot_mean_value/2); 
ascope_ylim_max = Vthreshold_absolute_envelope * 5; 
ylim(h_ax_ascope, [0 ascope_ylim_max]);


%% --- Main Simulation Loop ---
% For N_TARGETS mobile targets
fprintf('\nStarting continuous simulation...\n');
plot_update_interval = 100; 
plot_pause_duration = 0.001;
pulse_counter = 0; 

% While loop for continuous operation until figure is closed
while ishandle(main_fig)
    pulse_counter = pulse_counter + 1;
    current_time = (pulse_counter - 1) * dt; 
    
    antenna_pointing_angle = mod(scan_rate_rad * current_time, 2*pi);
    
    % Variables to store A-Scope data for the current pulse
    ascope_peak_amplitude = [];
    ascope_peak_apparent_range_km = [];

    % Loop over each target
    for k_target = 1:N_TARGETS
        % Move target k
        target_pos(k_target, :) = target_pos(k_target, :) + velocity(k_target, :) * dt; 
        
        if pulse_counter <= n_steps_alloc 
             target_pos_history{k_target}(pulse_counter, :) = target_pos(k_target, :);
        else 
            target_pos_history{k_target} = [target_pos_history{k_target}(2:end,:); target_pos(k_target,:)];
        end

        R_target_real_k = norm(target_pos(k_target, :)); 
        angle_target_real_k = atan2(target_pos(k_target, 2), target_pos(k_target, 1));
        
        % Reset target if it moves too far away
        if R_target_real_k > max_initial_range * 1.5 
            fprintf('Target %d reset — moved too far (R=%.1f km)...\n', k_target, R_target_real_k/1000);
            start_angle_rad_k = rand * 2*pi;
            start_range_k = (Rmax_na * 0.8) + rand * (max_initial_range - (Rmax_na * 0.8)); 
            if start_range_k > Rmax_na; speed_k = 600 + rand*300; else; speed_k = 400 + rand*250; end
            target_pos(k_target, :) = [start_range_k*cos(start_angle_rad_k), start_range_k*sin(start_angle_rad_k)];
            % Reset target direction always toward radar
            if norm(target_pos(k_target, :)) > 1e-3
                velocity_direction_k = -target_pos(k_target, :) / norm(target_pos(k_target, :));
            else
                random_direction_angle_k = rand * 2*pi;
                velocity_direction_k = [cos(random_direction_angle_k), sin(random_direction_angle_k)];
            end
            velocity(k_target, :) = speed_k * velocity_direction_k;
            
            target_pos_history{k_target}(:,:) = NaN; 
            if pulse_counter <= n_steps_alloc
                target_pos_history{k_target}(pulse_counter, :) = target_pos(k_target, :);
            else 
                 target_pos_history{k_target}(end, :) = target_pos(k_target, :); 
            end
            fprintf('New start Target %d: R=%.1f km, Ang=%.1f deg, Vel=%.1f m/s, Dir Toward Radar\n', k_target, start_range_k/1000, rad2deg(start_angle_rad_k), speed_k);
            continue; 
        end

        % Incorporating fluctuating targets (Swerling I/II)
        RCS_instantaneous_k = exprnd(RCS_mean_calc); 

        % Calculate received power and received voltage
        Prx_current_antenna_k = 0; S_current_receiver_k = 0; is_target_illuminated_k = false;
        V_s_pure_envelope_k = 0; 

        angular_difference_beam_k = abs(wrapToPi(antenna_pointing_angle - angle_target_real_k));
        if angular_difference_beam_k < (theta_bw_rad / 2) && R_target_real_k > R_min
            is_target_illuminated_k = true;
            Prx_num_calc_k = Pdg * Gant^2 * lambda^2 * RCS_instantaneous_k;
            Prx_den_calc_k = (4*pi)^3 * R_target_real_k^4;
            Prx_current_antenna_k = Prx_num_calc_k / Prx_den_calc_k;
            S_current_receiver_k = Prx_current_antenna_k / Ls;
            V_s_pure_envelope_k = sqrt(S_current_receiver_k); 
        end
        
        % Calculate exact round-trip signal travel time
        t_round_trip_k = 2 * R_target_real_k / c;

        % Add random noise and sum with signal
        sigma_noise_component_k = sqrt(Ntot_mean_value / 2);
        v_n_I_k = randn() * sigma_noise_component_k; v_n_Q_k = randn() * sigma_noise_component_k;
        random_signal_phase_k = 2*pi*rand;
        v_total_I_k = V_s_pure_envelope_k * cos(random_signal_phase_k) + v_n_I_k;
        v_total_Q_k = V_s_pure_envelope_k * sin(random_signal_phase_k) + v_n_Q_k;
        V_total_instantaneous_envelope_k = sqrt(v_total_I_k^2 + v_total_Q_k^2);

        % Cell-by-cell comparison against Vthreshold.
        % Targets in the blind zone do not appear on radar.
        is_in_physical_blind_zone_k = (R_target_real_k >= Rmax && R_target_real_k < Rmax_na);
        
        can_attempt_detection_k = is_target_illuminated_k && ...
                                ~is_in_physical_blind_zone_k && ...
                                R_target_real_k > R_min;
        
        if can_attempt_detection_k
            ascope_peak_amplitude(end+1) = V_total_instantaneous_envelope_k;
            ascope_peak_apparent_range_km(end+1) = mod(R_target_real_k, Rmax_na)/1000;

            if V_total_instantaneous_envelope_k > Vthreshold_absolute_envelope
                if t_round_trip_k < PRI
                    R_apparent_ppi_k = mod(R_target_real_k, Rmax_na);
                    detections(end+1, :) = [current_time, R_apparent_ppi_k, antenna_pointing_angle];
                end
            end
        end
    end % End of target loop

    % --- Plot Updates ---
    if mod(pulse_counter, plot_update_interval)==0
        % Update Real Position Display
        for k_target = 1:N_TARGETS
            valid_hist_idx = ~isnan(target_pos_history{k_target}(:,1));
            if any(valid_hist_idx)
                set(h_target_real_lines(k_target),'XData',target_pos_history{k_target}(valid_hist_idx,1),'YData',target_pos_history{k_target}(valid_hist_idx,2));
            else
                set(h_target_real_lines(k_target),'XData',NaN,'YData',NaN); 
            end
            set(h_target_curr_pos(k_target),'XData',target_pos(k_target,1),'YData',target_pos(k_target,2));
        end
        
        % Update PPI
        set(h_scan_line_ppi,'ThetaData',[antenna_pointing_angle antenna_pointing_angle]);
        if ~isempty(detections)
            scan_duration=360/scan_rate_deg; display_time_window=scan_duration*1.5; 
            valid_detections_to_plot=detections(detections(:,1)>=(current_time-display_time_window),:);
            if ~isempty(valid_detections_to_plot)
                set(h_detections_ppi_plot,'ThetaData',valid_detections_to_plot(:,3),'RData',valid_detections_to_plot(:,2));
                max_r_detected=max(valid_detections_to_plot(:,2));
                current_ppi_r_limit = ax_ppi.RLim(2);
                needed_ppi_r_limit=max([display_range_limit_initial,max_r_detected*1.05, Rmax_na*1.05]);
                 if current_ppi_r_limit < needed_ppi_r_limit || (max_r_detected < current_ppi_r_limit * 0.8 && current_ppi_r_limit > display_range_limit_initial)
                    ax_ppi.RLim(2)=needed_ppi_r_limit;
                    ppi_r_axis.TickValues=linspace(0,needed_ppi_r_limit,5);
                    ppi_r_axis.TickLabels=string(round(ppi_r_axis.TickValues/1000));
                 end
            else; set(h_detections_ppi_plot,'ThetaData',NaN,'RData',NaN); end
        else; set(h_detections_ppi_plot,'ThetaData',NaN,'RData',NaN); end
        
        current_max_real_dim_all_targets = 0;
        for k_target=1:N_TARGETS
             valid_hist_idx = ~isnan(target_pos_history{k_target}(:,1));
             if any(valid_hist_idx)
                max_coord_current_target_hist = max(abs(target_pos_history{k_target}(valid_hist_idx,:)), [], 'all');
                current_max_real_dim_all_targets = max(current_max_real_dim_all_targets, max_coord_current_target_hist);
             end
        end
        if current_max_real_dim_all_targets > 0
            current_xlim_real_val = get(h_ax_real, 'XLim');
            if current_max_real_dim_all_targets > current_xlim_real_val(2) || current_max_real_dim_all_targets < display_range_limit_initial 
                axis(h_ax_real, max(current_max_real_dim_all_targets * 1.1, display_range_limit_initial) * [-1 1 -1 1]);
            end
        else 
             axis(h_ax_real,display_range_limit_initial*[-1 1 -1 1]); 
        end

        % Update A-SCOPE Display (with multiple peaks if applicable)
        ascope_y_data = ones(1, num_ascope_points) * ascope_base_noise_amplitude .* (1 + 0.2*rand(1, num_ascope_points)); 
        if ~isempty(ascope_peak_apparent_range_km)
            for p_idx = 1:length(ascope_peak_apparent_range_km)
                r_peak_km = ascope_peak_apparent_range_km(p_idx);
                amp_peak = ascope_peak_amplitude(p_idx);
                if r_peak_km <= ascope_range_display_km 
                    [~, idx_target_ascope] = min(abs(ascope_range_bins_km - r_peak_km));
                    pulse_width_display_ascope = max(1, round( (R_res/1000) / (ascope_range_display_km/num_ascope_points) )); 
                    idx_start = max(1, idx_target_ascope - floor(pulse_width_display_ascope/2));
                    idx_end = min(num_ascope_points, idx_target_ascope + ceil(pulse_width_display_ascope/2)-1);
                    
                    if idx_start <= idx_end 
                        num_pts_rise = ceil((idx_end-idx_start+1)/2);
                        if num_pts_rise > 0
                            pulse_shape_rise = linspace(ascope_y_data(idx_start), amp_peak, num_pts_rise);
                            ascope_y_data(idx_start : idx_start+num_pts_rise-1) = max(ascope_y_data(idx_start : idx_start+num_pts_rise-1), pulse_shape_rise);
                            
                            num_pts_fall = (idx_end-idx_start+1) - num_pts_rise;
                            if num_pts_fall > 0
                                pulse_shape_fall = linspace(amp_peak, ascope_y_data(idx_end), num_pts_fall+1);
                                ascope_y_data(idx_start+num_pts_rise : idx_end) = max(ascope_y_data(idx_start+num_pts_rise : idx_end), pulse_shape_fall(2:end));
                            end
                        end
                    end
                end
            end
        end
        set(h_ascope_plot, 'XData', ascope_range_bins_km, 'YData', ascope_y_data);
        
        drawnow limitrate; 
        pause(plot_pause_duration);
    end

    % Simulation Stop Condition (by time limit or closing the window)
    if current_time >= sim_time && sim_time > 0 
        fprintf('Simulation time (%.1f s) reached. Terminating.\n', sim_time);
        break; 
    end
end % End of Main Loop (while)

if ishandle(main_fig)
    fprintf('Simulation stopped by user or time limit reached.\n');
else
    fprintf('Simulation window closed by user.\n');
end