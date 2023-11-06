function [step_stats, angles, moments, grf, startStance] = ProcessLab3Data(filename)
    % Load and Process Data
    loadedData = load(filename);
    freq = 100;
    forceplate = downsample(loadedData.ForcePlate_Data(:,1:6), 10);
    COP = downsample(loadedData.ForcePlate_Data(:,7:9),10); 

    % Filtering data
    [a,b] = butter(2,10/50);
    forceplate = filtfilt(a,b,forceplate);
    COP = filtfilt(a,b,COP);
    [a,b] = butter(2,5/50);
    markerpos = filtfilt(a,b,loadedData.Marker_Data);

    % Define ground reaction force with our reference frame (opposite z)
    grf = forceplate(:,1:2);
    grf(:,3) = -forceplate(:,3);

    % Calculate threshold for start of stance phase based on grf
    noForce = grf(1:35,3);
    mean_noForce = mean(noForce);
    std_grf = std(noForce);
    startStance = mean_noForce + 4*std_grf;

    COP = -COP/1000; % Convert from mm to m and rotate reference frame 180 degrees

    % marker positions in meters (separate [x,y,z] coords of each marker)
    % Then,reverse x and y direction to match reference frame
    Rankle = markerpos(:,1:3)/1000;
    Rankle(:,1:2) = -Rankle(:,1:2);
    Rtip = markerpos(:,4:6)/1000;
    Rtip(:,1:2) =-Rtip(:,1:2);
    Rheel = markerpos(:,7:9)/1000;
    Rheel(:,1:2) = -Rheel(:,1:2);
    Rknee = markerpos(:,10:12)/1000;
    Rknee(:,1:2) =-Rknee(:,1:2);
    Rtroch = markerpos(:,13:15)/1000;
    Rtroch(:,1:2) =-Rtroch(:,1:2);
    Lankle = markerpos(:,16:18)/1000;
    Lankle(:,1:2) =-Lankle(:,1:2);
    Ltip = markerpos(:,19:21)/1000;
    Ltip(:,1:2) =-Ltip(:,1:2);
    Lheel = markerpos(:,22:24)/1000;
    Lheel(:,1:2) =-Lheel(:,1:2);
    Lknee = markerpos(:,25:27)/1000;
    Lknee(:,1:2) =-Lknee(:,1:2);
    Ltroch = markerpos(:,28:30)/1000;
    Ltroch(:,1:2) =-Ltroch(:,1:2);
    Rshldr = markerpos(:,31:33)/1000;
    Rshldr(:,1:2) =-Rshldr(:,1:2);
    Lshldr = markerpos(:,34:36)/1000;
    Lshldr(:,1:2) =-Lshldr(:,1:2);
    Rback = markerpos(:,37:39)/1000;
    Rback(:,1:2) =-Rback(:,1:2);
    clear('markerpos')
    

    % Define segment angles for right leg points
    Rtrunk_angles = getAngles(Rtroch, Rshldr);
    Rthigh_angles = getAngles(Rknee, Rtroch);
    Rshank_angles = getAngles(Rankle, Rknee);
    Rfoot_angles = pi + getAngles(Rheel, Rtip);

    % Take second derivatives of theta wrt time rad/s^2
    dt = 1/freq;
    aa_Rthigh = Deriv2(Rthigh_angles, dt, 1); 
    aa_Rshank = Deriv2(Rshank_angles, dt, 1); 
    aa_Rfoot = Deriv2(Rfoot_angles, dt, 1); 

    % Dictionary containing key anthropomorphic data
    Anth = struct('thigh_com', 0.3612, 'shank_com', 0.4416, 'foot_com', 0.4014, ...
                  'thigh_m', 0.1478, 'shank_m', 0.0481, 'foot_m', 0.0129, ...
                  'thigh_r', 0.369, 'shank_r', 0.271, 'foot_r', 0.299);
    mBody = 63; % Subject total body mass [kg]

    % Calculate COM for each segment (Right leg) for each frame
    [rthigh_COM, rthigh_len] = GetSegCOMs(Rtroch, Rknee, Rthigh_angles, Anth.thigh_com);
    [rshank_COM, rshank_len] = GetSegCOMs(Rknee, Rankle, Rshank_angles, Anth.shank_com);
    [rfoot_COM, rfoot_len] = GetSegCOMs(Rheel, Rtip, Rfoot_angles, Anth.foot_com);

    % Calculate mass for each segment
    mThigh = mBody * Anth.thigh_m;
    mShank = mBody * Anth.shank_m;
    mFoot = mBody * Anth.foot_m;

    % Calculate radius of gyration for each segment
    radThigh = rthigh_len * Anth.thigh_r;
    radShank = rshank_len * Anth.shank_r;
    radFoot = rfoot_len * Anth.foot_r;

    % Calculate I for each segment
    Ithigh = mThigh * radThigh^2;
    Ishank = mShank * radShank^2;
    Ifoot = mFoot * radFoot^2;

    % Calculate Linear Accel for each segment COM
    la_rthigh = [Deriv2(rthigh_COM, dt, 1), Deriv2(rthigh_COM, dt, 2)];
    la_rshank = [Deriv2(rshank_COM, dt, 1), Deriv2(rshank_COM, dt, 2)];
    la_rfoot = [Deriv2(rfoot_COM, dt, 1), Deriv2(rfoot_COM, dt, 2)];

    % Use Inverse Dynamics Functions to Calculate Reaction Forces and Moments
    [Ra, Ma] = footdynamics(mFoot, la_rfoot, grf, COP, rfoot_COM, aa_Rfoot, Ifoot, Rankle);
    [Rk, Mk] = shankdynamics(mShank, la_rshank, rshank_COM, aa_Rshank, Ishank, Rankle, Rknee, Ra, Ma);
    [~, Mh] = thighnamics(mThigh, la_rthigh, rthigh_COM, aa_Rthigh, Ithigh, Rknee, Rtroch, Rk, Mk);

    % Calculate gait speed, step length, cadence
    % Let's use average right knee x-velocity for gait speed
    v = zeros(size(Rknee, 1) - 1, 1);
    for j = 1:size(Rknee, 1) - 1
        v(j) = (Rknee(j + 1, 1) - Rknee(j, 1)) / dt;
    end
    gait_speed = mean(v); % m/s

    % Find step length by finding left and right heel contact
    [~, rhc_ind] = FindFirstMin(Rheel(:,3));
    [~, lhc_ind] = FindFirstMin(Lheel(:,3));
    step_len = Lheel(lhc_ind, 1) - Rheel(rhc_ind, 1); % meters

    cadence = 60 * (gait_speed / step_len); % steps per min

    % Check heel distance results visually
    % hold on
    % plot(grf(:, 3) / 500)
    % plot(Rheel(:, 3))
    % plot(Lheel(:, 3))
    % xline(rhc_ind)
    % xline(lhc_ind)
    % legend('grf', 'right heel', 'left heel')

    % Results to return
    step_stats = [gait_speed, step_len, cadence];
    angles = [Rfoot_angles, Rshank_angles, Rthigh_angles];
    moments = [Ma, Mk, Mh];
    % Also return grf and startStance

    % Sub-functions to keep everything nice and clean!
    function theta = getAngles(L, U)
        % L and U should be nx3 arrays with [x, y, z] coordinates of joint markers.
        % The function will then iterate through the points and calculate an angle
        % vector which is (nx1) for each angle in the gait cycle. Function uses
        % atan2d to return theta as an array of angles in degrees where theta is
        % between -90 and 90. L is the lower point, U is the upper point.
        if size(L) ~= size(U)
            error('Input arrays must be the same size');
        end

        theta = zeros(size(L, 1), 1);

        for i = 1:size(L, 1)
            theta(i) = atan2((U(i, 3) - L(i, 3)), (U(i, 1) - L(i, 1)));
        end

    end

    function [a] = Deriv2(pos, dt, dir)
        % Takes Position input, which should be an nx1:3 position array and dir
        % which defines the direction derivative is to be taken wrt. dt is a
        % constant change in time defined based on the sample rate.

        a = zeros(size(pos, 1), 1);
        x = pos(:, dir);

        for i = 2:size(pos, 1) - 1
            a(i) = (x(i + 1) - 2 * x(i) + x(i - 1)) / (dt^2);
        end

    end

    function [seg_COM, seg_len] = GetSegCOMs(U, L, angles, d)
        % This function takes position data for upper and lower join and uses
        % angles to calculate average segment length as well as the COM using
        % anthropomorphism constants, d.

        seg_len = sqrt((U(:, 1) - L(:, 1)).^2 + (U(:, 3) - L(:, 3)).^2); % calculate the segment length in each frame
        seg_len = mean(seg_len); % assign the average value as a constant for the segment length

        x_com = U(:, 1) - d * seg_len * cos(angles);
        z_com = U(:, 3) - d * seg_len * cos(angles);

        seg_COM = [x_com, z_com];

    end

    % Inverse Dynamics functions written for ME21 with minor adjustments to
    % match the dimensions of the data being used for this lab:
    function [Ra, Ma] = footdynamics(m, a, grf, cop, com, alpha, I, ankle)
        % Inputs: m = segment mass scalar
        % a = linear acceleration of COM [ax, az]
        % grf = ground reaction force matrix [gx, gy, gz]
        % cop = center of pressure, nx3 = [x, y, z]
        % com = center of mass nx2 = [x, z]
        % alpha = angular acceleration of segment nx1
        % I = moment of inertia of segment, scalar
        % ankle = ankle position nx3 = [x, y, z]

        Ra(:, 1) = m * a(:, 1) - grf(:, 1);
        Ra(:, 2) = m * a(:, 2) + 9.81 * m - grf(:, 3);

        M1 = grf(:, 3) .* (cop(:, 1) - com(:, 1));
        M2 = grf(:, 1) .* com(:, 2);
        M3 = Ra(:, 2) .* (com(:, 1) - ankle(:, 1));
        M4 = Ra(:, 1) .* (ankle(:, 3) - com(:, 2));

        Ma = I * alpha - M1 - M2 + M3 + M4;
    end

    function [Rk, Mk] = shankdynamics(m, a, com, alpha, I, ankle, knee, Ra, Ma)
        % Inputs for shank, but same as for foot, except now with Ra, Ma as
        % inputs. Be careful about nx2 and nx3 vectors.

        Rk(:, 1) = m * a(:, 1) + Ra(:, 1);
        Rk(:, 2) = m * a(:, 2) + Ra(:, 2) + 9.81 * m;

        M1 = Rk(:, 2) .* (com(:, 1) - knee(:, 1));
        M2 = Rk(:, 1) .* (knee(:, 3) - com(:, 2));
        M3 = Ra(:, 2) .* (ankle(:, 1) - com(:, 1));
        M4 = Ra(:, 1) .* (com(:, 2) - ankle(:, 3));

        Mk = I * alpha + Ma + M1 + M2 + M3 + M4;
    end

    function [Rh, Mh] = thighnamics(m, a, com, alpha, I, knee, hip, Rk, Mk)
        % Outputs: hip reaction and moment vector
        Rh(:, 1) = m * a(:, 1) + Rk(:, 1);
        Rh(:, 2) = m * a(:, 2) + Rk(:, 2) + 9.81 * m;

        M1 = Rh(:, 2) .* (com(:, 1) - hip(:, 1));
        M2 = Rh(:, 1) .* (hip(:, 3) - com(:, 2));
        M3 = Rk(:, 2) .* (knee(:, 1) - com(:, 1));
        M4 = Rk(:, 1) .* (com(:, 2) - knee(:, 3));

        Mh = I * alpha + Mk + M1 + M2 + M3 + M4;
    end

    function [min_value, min_index] = FindFirstMin(y)
        % Initialize minimum value and index
        min_value = inf;
        min_index = -1;

        % Iterate through the vector y
        for i = 2:(length(y) - 1)
            if y(i) < y(i - 1) && y(i) < y(i + 1)
                % Found a local minimum with larger values on both sides
                if y(i) < min_value
                    min_value = y(i);
                    min_index = i;
                    return
                end
            end
        end

        % If a valid minimum was found, min_index will be greater than 0
        if min_index == -1
            error('No valid local minimum found in the input vector.');
        end
    end
end
