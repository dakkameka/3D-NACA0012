% THE PLOTTING CODE TO RUN AFTER TERMINATING SOMETHING EARLY

close all;
plot_variable = 'all';

% plotStreamlines(cell_centers, Q_numerical, face_centers_x, face_centers_y, BC_airfoil_faces)

%% SHEAR STRESS AND SKIN FRICTION COEFFICIENT ALONG AIRFOIL
%% WALL SHEAR STRESS FROM CELL-CENTERED GRADIENTS

% if DIRECT_NS
%     disp("plotting wall shear stress and skin friction c_f from cell center gradients...")
% 
%     mu_w = mu_inf;
%     airfoil_faces = BC_airfoil_faces;
%     x_fc = face_centers(airfoil_faces,1);
%     y_fc = face_centers(airfoil_faces,2);
%     n_faces = length(airfoil_faces);
% 
%     tau_w_cell = zeros(n_faces,1);
%     cf_cell = zeros(n_faces,1);
% 
%     for k = 1:n_faces
%         f = airfoil_faces(k);
% 
%         % OUTWARD NORMAL AT WALL FACE
%         n_vec = [nx(f), ny(f)];
%         t_vec = [-ny(f), nx(f)];
% 
%         % WALL CELL (left of wall face)
%         c = cell_1_of_face(f);
% 
%         % CELL-CENTER GRADIENTS OF Q
%         dQdx_cell = dQdx(:,c);
%         dQdy_cell = dQdy(:,c);
% 
%         % CELL VELOCITY
%         Qc = Q_numerical(:,c);
%         rho_c = Qc(1);
%         u_c = Qc(2)/rho_c;
%         v_c = Qc(3)/rho_c;
% 
%         % VELOCITY GRADIENTS AT CELL CENTER
%         du_dx = (dQdx_cell(2) - u_c * dQdx_cell(1)) / rho_c;
%         du_dy = (dQdy_cell(2) - u_c * dQdy_cell(1)) / rho_c;
%         dv_dx = (dQdx_cell(3) - v_c * dQdx_cell(1)) / rho_c;
%         dv_dy = (dQdy_cell(3) - v_c * dQdy_cell(1)) / rho_c;
% 
%         n_hat = n_vec / norm(n_vec);
%         t_hat = t_vec / norm(t_vec);
% 
%         GradU = [du_dx, du_dy; dv_dx, dv_dy];
%         dudn_vec = GradU * n_hat';
%         dutdn = dot(t_hat, dudn_vec);
% 
%         tau_w_cell(k) = abs(mu_w * dutdn);
%         cf_cell(k) = tau_w_cell(k) / (0.5 * rho_inf * V_inf^2);
%     end
% 
%     % SPLIT TOP/BOTTOM AND PLOT AS BEFORE
%     x_le = min(x_fc); x_te = max(x_fc); chord_length = x_te - x_le;
%     x_c = (x_fc - x_le) / chord_length;
%     is_top = y_fc >= 0; is_bot = ~is_top;
%     [xtop, itop] = sort(x_c(is_top)); cf_top = cf_cell(is_top); cf_top = cf_top(itop);
%     [xbot, ibot] = sort(x_c(is_bot)); cf_bot = cf_cell(is_bot); cf_bot = cf_bot(ibot);
% 
%     figure; hold on;
%     plot(xtop, cf_top, 'b-', 'LineWidth', 2);
%     plot(xbot, cf_bot, 'r-', 'LineWidth', 2);
%     xlabel('x/c'); ylabel('c_f (cell grad)');
%     title('Skin Friction Coefficient');
%     legend('Top Surface', 'Bottom Surface');
%     grid on;
% end
% 


    % plotlimit = 20; figure; hold on;
    % plot([x1'; x2'], [y1'; y2'], 'k-');  % faces
    % if plotNorms == true
    %     quiver(face_centers_x, face_centers_y, nx, ny, 0.3, 'r');  % face normals
    % end
    % scatter(face_centers_x, face_centers_y, 5, 'b', 'filled');   % face centers
    % scatter(cell_centers_x, cell_centers_y, 10, 'g', 'filled');  % cell centers
    % scatter(ghost_centers(:,1), ghost_centers(:,2), 15, 'm', 'filled'); % ghost cell centers in magenta
    % 
    % xlim([-plotlimit plotlimit]); ylim([-plotlimit plotlimit]);
    % xlabel('x'); ylabel('y');
    % title('Face Normals, Face Centers, Cell Centers, and Ghost Cells');
    % legend('Faces', 'Normals', 'Face Centers', 'Cell Centers', 'Ghost Cell Centers');



    disp("Moved on to plotting may take awhile heehee")

    rho_vals = Q_numerical(1,:);
    u_vals = Q_numerical(2,:) ./ rho_vals;   v_vals = Q_numerical(3,:) ./ rho_vals;
    E_vals = Q_numerical(4,:);    p_vals = (gamma - 1) * (E_vals - 0.5 * rho_vals .* (u_vals.^2 + v_vals.^2));
    T_vals = p_vals ./ (rho_vals * 286.69);    velocity_mag = sqrt(u_vals.^2 + v_vals.^2);
    mach_vals = velocity_mag ./ sqrt(gamma * p_vals ./ rho_vals); 

    variables = {'density', 'u', 'v', 'pressure', 'temperature', 'velocity', 'mach'};
    titles = {'Density Distribution', 'X-Velocity Distribution', 'Y-Velocity Distribution',...
          'Pressure Distribution', 'Temperature Distribution', 'Velocity Magnitude', 'Mach Number'};
    data = {rho_vals, u_vals, v_vals, p_vals, T_vals, velocity_mag, mach_vals};

    if strcmpi(plot_variable, 'all')
        plot_indices = 1:length(variables);
    else
        plot_indices = find(strcmpi(variables, plot_variable));
    end

    for idx = plot_indices
        figure('Name', titles{idx}); hold on; title(titles{idx});
        xlabel('x'); ylabel('y'); colormap jet; colorbar; axis equal;

        for i = 1:n_cells
            faces_of_cell = faces_per_cell{i};
            xv = [x1(faces_of_cell); x2(faces_of_cell)];
            yv = [y1(faces_of_cell); y2(faces_of_cell)];
            k = convhull(xv, yv);
            fill(xv(k), yv(k), data{idx}(i), 'EdgeColor', 'none');
        end

        % Plot mesh edges
        %plot([x1'; x2'], [y1'; y2'], 'k-', 'LineWidth', 0.1);
        hold off;
        xlim([-plotlimit plotlimit]); ylim([-plotlimit plotlimit]);
    end


figure('Name', 'Residual Convergence');semilogy(residual_history, 'LineWidth', 2);
title('Residual Convergence History');xlabel('Iteration');
ylabel('Residual (log scale)');legend('L2 Norm Residual');


        disp("the mach number and cp on airfoil surface now :D");

        V_inf = sqrt(u_inf^2 + v_inf^2);

        % extract airfoil face centers that we already found 
        airfoil_faces = BC_airfoil_faces;
        x_fc = face_centers(airfoil_faces,1); y_fc = face_centers(airfoil_faces,2);
        n_faces = length(airfoil_faces);

        % normalizing and mapping it to 0 to 1 
        x_le = min(x_fc); x_te = max(x_fc); 
        chord_length = x_te - x_le;x_c = (x_fc - x_le) / chord_length;

        cp_vals = zeros(n_faces,1);
        mach_vals = zeros(n_faces,1);

        for i = 1:n_faces
            f = airfoil_faces(i);
            c = cell_1_of_face(f); % interior cell
            Q = Q_numerical(:, c);
            rho = Q(1);
            u = Q(2)/rho; v = Q(3)/rho;
            E = Q(4);
            p = (gamma - 1)*(E - 0.5*rho*(u^2 + v^2));
            a = sqrt(gamma * p / rho);
            V = sqrt(u^2 + v^2);
            cp_vals(i) = (p - p_inf) / (0.5 * rho_inf * V_inf^2);
            mach_vals(i) = V / a;
        end

        % SPLIT TOP AND BOTTOM BY Y
        is_top = y_fc >= 0;
        is_bot = ~is_top;

        % SORT TOP
        [xtop, itop] = sort(x_c(is_top));
        cp_top = cp_vals(is_top); cp_top = cp_top(itop);
        mach_top = mach_vals(is_top); mach_top = mach_top(itop);

        % SORT BOTTOM
        [xbot, ibot] = sort(x_c(is_bot));
        cp_bot = cp_vals(is_bot); cp_bot = cp_bot(ibot);
        mach_bot = mach_vals(is_bot); mach_bot = mach_bot(ibot);

        % PLOT Cp
        figure; hold on;
        plot(xtop, cp_top, 'b-', 'LineWidth', 2); plot(xbot, cp_bot, 'r-', 'LineWidth', 2);
        set(gca, 'YDir','reverse'); xlabel('x/c'); ylabel('Cp'); title('Pressure Coefficient on Airfoil');
        legend('Top Surface', 'Bottom Surface'); 

        % PLOT Mach
        figure; hold on;
        plot(xtop, mach_top, 'b-', 'LineWidth', 2); plot(xbot, mach_bot, 'r-', 'LineWidth', 2);
        xlabel('x/c'); ylabel('Mach Number'); title('Mach Number on Airfoil');
        legend('Top Surface', 'Bottom Surface'); grid on;


    disp("smooth contour plots heehee theyre more aesthetic :D");

    % COMPUTE FLOW VARIABLES
    rho_vals = Q_numerical(1,:);
    u_vals = Q_numerical(2,:) ./ rho_vals;
    v_vals = Q_numerical(3,:) ./ rho_vals;
    E_vals = Q_numerical(4,:);
    p_vals = (gamma - 1) * (E_vals - 0.5 * rho_vals .* (u_vals.^2 + v_vals.^2));
    velocity_mag = sqrt(u_vals.^2 + v_vals.^2);
    mach_vals = velocity_mag ./ sqrt(gamma * p_vals ./ rho_vals);

    % TRIANGULATION BASED ON CELL CENTERS
    tri = delaunay(cell_centers_x, cell_centers_y);
    tri_centers_x = mean(cell_centers_x(tri), 2);
    tri_centers_y = mean(cell_centers_y(tri), 2);

    % AIRFOIL POLYGON FROM FACE CENTERS
    airfoil_x = face_centers(BC_airfoil_faces,1);
    airfoil_y = face_centers(BC_airfoil_faces,2);
    k = convhull(airfoil_x, airfoil_y);
    airfoil_boundary_x = airfoil_x(k);
    airfoil_boundary_y = airfoil_y(k);

    % REMOVE TRIANGLES INSIDE AIRFOIL
    inside = inpolygon(tri_centers_x, tri_centers_y, ...
                       airfoil_boundary_x, airfoil_boundary_y);
    tri = tri(~inside, :);

    % DATA TO PLOT
    contour_vars = {rho_vals, p_vals, velocity_mag, mach_vals};
    contour_titles = {'Density Contour', 'Pressure Contour', ...
                      'Velocity Magnitude Contour', 'Mach Number Contour'};

    % LOOP THROUGH SELECTED VARIABLES
    for i = 1:length(contour_vars)
        figure('Name', contour_titles{i});
        trisurf(tri, cell_centers_x, cell_centers_y, contour_vars{i}, ...
                'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2); axis equal tight; colormap jet; colorbar;
        title(contour_titles{i});
        xlabel('x'); ylabel('y');
        xlim([-plotlimit plotlimit]);
        ylim([-plotlimit plotlimit]); grid off;
    end



    disp("plotting quiver plot :D")

    % COMPUTE FLOW VARIABLES
    rho_vals = Q_numerical(1,:);
    u_vals = Q_numerical(2,:) ./ rho_vals;
    v_vals = Q_numerical(3,:) ./ rho_vals;
    velocity_mag = sqrt(u_vals.^2 + v_vals.^2);

    % GENERATE GRID FOR QUIVERS
    xg = linspace(min(cell_centers_x), max(cell_centers_x), 1000);
    yg = linspace(min(cell_centers_y), max(cell_centers_y), 3500);
    [X, Y] = meshgrid(xg, yg);
    U = griddata(cell_centers_x, cell_centers_y, u_vals, X, Y, 'linear');
    V = griddata(cell_centers_x, cell_centers_y, v_vals, X, Y, 'linear');

    % MASK OUT POINTS INSIDE AIRFOIL
    airfoil_x = face_centers(BC_airfoil_faces,1);
    airfoil_y = face_centers(BC_airfoil_faces,2);
    k = convhull(airfoil_x, airfoil_y);
    mask = inpolygon(X, Y, airfoil_x(k), airfoil_y(k));
    U(mask) = NaN;
    V(mask) = NaN;

    figure('Name', 'Velocity Magnitude with Quivers'); hold on;
    xlabel('x'); ylabel('y'); title('Velocity Magnitude and Velocity Vectors');
    colormap jet; colorbar; axis equal;

    % PLOT VELOCITY MAGNITUDE
    for i = 1:n_cells
        faces_of_cell = faces_per_cell{i};
        xv = [x1(faces_of_cell); x2(faces_of_cell)];
        yv = [y1(faces_of_cell); y2(faces_of_cell)];
        k2 = convhull(xv, yv);
        fill(xv(k2), yv(k2), velocity_mag(i), 'EdgeColor', 'none');
    end

    % PLOT QUIVERS SMALL AND BLACK
    quiver(X, Y, U, V, 0.99, 'k'); % 0.5 = scale down arrow length

    xlim([-plotlimit plotlimit]); ylim([-plotlimit plotlimit]);


    function plotStreamlines(cell_centers, Q, face_centers_x, face_centers_y, BC_airfoil_faces)

    % GET VELOCITY COMPONENTS (u,v) FROM CONSERVED VARS
    rho = Q(1,:)';
    u   = Q(2,:)' ./ rho;
    v   = Q(3,:)' ./ rho;

    % GRID FOR INTERPOLATION
    x_cell = cell_centers(:,1);
    y_cell = cell_centers(:,2);
    xq = linspace(min(x_cell), max(x_cell), 600);
    yq = linspace(min(y_cell), max(y_cell), 600);
    [Xq, Yq] = meshgrid(xq, yq);

    % INTERPOLATE VELOCITIES TO STRUCTURED GRID
    Uq = griddata(x_cell, y_cell, u, Xq, Yq, 'natural');
    Vq = griddata(x_cell, y_cell, v, Xq, Yq, 'natural');

    % REMOVE ANY NANs FROM INTERPOLATION
    Uq(isnan(Uq)) = 0;
    Vq(isnan(Vq)) = 0;

    % MASK VELOCITIES INSIDE AIRFOIL SO STREAMLINES AVOID IT
    airfoil_x = face_centers_x(BC_airfoil_faces);
    airfoil_y = face_centers_y(BC_airfoil_faces);
    k = convhull(airfoil_x, airfoil_y);
    airfoil_x = airfoil_x(k);
    airfoil_y = airfoil_y(k);
    mask = inpolygon(Xq, Yq, airfoil_x, airfoil_y);
    Uq(mask) = NaN;
    Vq(mask) = NaN;

    % DEFINE STARTING STREAMLINE POINTS
    % % % % % x_start = linspace(min(xq), min(xq) + 0.1*(max(xq)-min(xq)), 180);
    y_start = linspace(min(yq), max(yq), 120);
    %y_start = linspace(-6, 2, 60);   
    [startX, startY] = meshgrid(x_start, y_start);

    % PLOT STREAMLINES
    figure;
    streamline(Xq, Yq, Uq, Vq, startX, startY)
    % hold on

    % FILL AIRFOIL SHAPE WITH WHITE AND OUTLINE IN BLACK
    fill(airfoil_x, airfoil_y, 'w', 'EdgeColor', 'k', 'LineWidth', 1.5);

    axis equal tight
    xlabel('x'); ylabel('y');
    title('STREAMLINES AROUND AIRFOIL');
    set(gcf, 'Color', 'w');

end