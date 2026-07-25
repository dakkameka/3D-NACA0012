% attempt at DIRECT NAVIER STOKES
  % and first order DNS

% speed this baby up



%% USER SETUP__________________________________________________________________________________________________________________________________
restart = false;

CFL = 1;
mach = 0.5;
order = 2;
alpha = 0;
timeStep = 1; % 1 for global 2 for local... USE GLOBAL FOR ORDER 2
orderDesired = 8; % like 6 orders, 7 orders, or 8 orders
plotMesh = false;
plotStreamlines = false;
plotFill = true;
plotContour = true;
plotNorms = false; 
plotCpMach = false;
filename = 'grid_file.in';
plotlimit = 2;
residual = inf;
%residual = 1*10^(-10);
tolerance = 0;
max_iterations = 1000000;
plot_variable = 'all'; %  'u', 'v', 'temperature', 'pressure', 'velocity', or 'all' 
limiter = true;

rho = 1.4; p = 1; gamma = 1.4; c = sqrt( gamma * p / rho); VMAG = mach * c; u = VMAG*cosd(alpha); v = VMAG*sind(alpha);
rho_inf = rho; p_inf   = p; v_inf   = v; c_inf = sqrt(gamma * p_inf / rho_inf); u_inf   = mach * c_inf; E_inf   = p_inf / (gamma - 1) + 0.5 * rho_inf * (u_inf^2 + v_inf^2); rho_inf = rho; p_inf   = p; v_inf   = v; c_inf = c; u_inf   = u;V_inf = sqrt(u_inf^2 + v_inf^2);

DIRECT_NS = true; %____________________________________________________
R = 286.69;
Pr = 0.72;
gamma = 1.4;
Cp = gamma*R / (gamma-1);
Re = 5000;
Lref = 1;
mu_inf = rho_inf * VMAG * Lref / Re;
k_inf = mu_inf * Cp / Pr;

%% INITIALIZING__________________________________________________________________________________________________________________________________
if restart == true % restarts :)
    iteration = 1; residual = inf;
    tolerance = 1*10^(-8); residual_history = [];
    CD_history = []; CL_history = [];
    [n_cells, ~, n_points, cell_centers, nx, ny, face_centers, ...
     x1, x2, y1, y2, x, y, lengths, cell_1_of_face, cell_2_of_face] = loadGridData('grid_file.in');
    faces_per_cell = cellsFaces(n_cells, cell_1_of_face, cell_2_of_face);
    face_centers_x = face_centers(:,1); face_centers_y = face_centers(:,2);
    cell_centers_x = cell_centers(:,1); cell_centers_y = cell_centers(:,2);
    Vol = computeCellVolumes(n_cells, faces_per_cell, x1, y1, x2, y2, cell_centers);
    [BC_inlet_faces, BC_inlet_cells] = getInletBoundary(cell_1_of_face, cell_2_of_face, face_centers);
    [BC_airfoil_faces, BC_airfoil_cells] = getAirfoilBoundary(cell_1_of_face, cell_2_of_face, face_centers);
    [BC_outlet_faces, BC_outlet_cells] = getOutletBoundary(cell_1_of_face, cell_2_of_face, face_centers);
    [BC_farfield_faces, BC_farfield_cells] = getFarfieldBoundary(cell_1_of_face, cell_2_of_face, BC_inlet_faces, BC_airfoil_faces, BC_outlet_faces);
    ghost_centers = computeGhostCenters(BC_airfoil_faces, cell_1_of_face, nx, ny, face_centers, cell_centers);
    Q_numerical = zeros(4, n_cells); % INITIALIZE ALL CELLS TO FREESTREAM
    Q_numerical(1,:) = rho_inf;                 Q_numerical(2,:) = rho_inf * u_inf;
    Q_numerical(3,:) = rho_inf * v_inf;         Q_numerical(4,:) = E_inf;
end

%% ITERATION LOOP___________________________________________________________________________________________________________________________________________________________
while residual > tolerance && iteration < max_iterations
    dt = computeDT(Q_numerical, CFL, Vol, faces_per_cell, nx, ny, lengths, gamma, n_cells, timeStep);

    % step 1 of SSP RK2
    [QL, QR, dQdx, dQdy] = reconstruct(Q_numerical, cell_1_of_face, cell_2_of_face, nx, ny, BC_inlet_faces, BC_airfoil_faces, BC_outlet_faces, BC_farfield_faces, mach, gamma, alpha, cell_centers, face_centers, faces_per_cell,order);
        if order == 2
        [dQdx_f, dQdy_f] = computeGradientsAtFaces(Q_numerical, dQdx, dQdy, cell_centers, face_centers, cell_1_of_face, cell_2_of_face, nx, ny, BC_airfoil_faces, 'adiabatic', []);  % or 'isothermal', 300
        Fv = computeViscousFluxes(Q_numerical, dQdx_f, dQdy_f, nx, ny, cell_1_of_face, cell_2_of_face, mu_inf, k_inf, gamma, R);
        end
    [Fluxes] = computeFluxes(QL, QR, nx, ny, gamma,dQdx_f, dQdy_f, Fv); 
    Res = computeResiduals(Fluxes, cell_1_of_face, cell_2_of_face, lengths, n_cells, Vol);
    Q_star = Q_numerical + dt .* (Res);

    % step 2 of SSP RK2
    [QL, QR, dQdx, dQdy] = reconstruct(Q_star, cell_1_of_face, cell_2_of_face, nx, ny, BC_inlet_faces, BC_airfoil_faces, BC_outlet_faces, BC_farfield_faces, mach, gamma, alpha, cell_centers, face_centers, faces_per_cell,order);
        if order ==2
        [dQdx_f, dQdy_f] = computeGradientsAtFaces(Q_star, dQdx, dQdy, cell_centers, face_centers, cell_1_of_face, cell_2_of_face, nx, ny, BC_airfoil_faces, 'adiabatic', []);  % or 'isothermal', 300
        Fv = computeViscousFluxes(Q_star, dQdx_f, dQdy_f, nx, ny, cell_1_of_face, cell_2_of_face, mu_inf, k_inf, gamma, R);
        end
    [Fluxes1] = computeFluxes(QL, QR, nx, ny, gamma,dQdx_f, dQdy_f, Fv);
    Res = computeResiduals(Fluxes1, cell_1_of_face, cell_2_of_face, lengths, n_cells, Vol);
    Q_new = 0.5 * (Q_numerical + Q_star + dt .* (Res));

    Q_numerical = Q_new;    iteration = iteration + 1;
    %Q_ghost = computeGhostStates(Q_numerical, BC_airfoil_faces, cell_1_of_face, face_centers, ghost_centers, nx, ny);
     
    Res_combined = sqrt(sum(Res.^2, 1));  residual = norm(Res_combined, 2);


        residual_history = [residual_history, residual]; 
    
        if mod(iteration,10) == 0
        if max_iterations ~= Inf % if finite number of itterantions keep length to max
            iter_format_str = sprintf("%d", floor(log10(max_iterations)));
        else % if inifity itterations let it dynamically change length
            iter_format_str = "";
        end
        print_format = sprintf("iteration: %%%sd | residual: %%.5e | tolerance: %%.2e \\n", iter_format_str);
        fprintf(print_format, iteration, residual, tolerance);
        %[CD, CL] = computeForceCoefficients(Q_numerical, BC_airfoil_faces, cell_1_of_face, nx, ny, lengths, gamma, rho_inf, V_inf, alpha, 1);
        %CD_history = [CD_history, CD];
        %CL_history = [CL_history, CL];

        %fprintf('Cd = %.6f | Cl = %.6f\n', CD, CL);
        end


    tolerance = residual_history(1) * 10^(-orderDesired);
end



%% PLOTTING_________________________________________________________________________________________________________________________________________________________________
        if plotMesh == true
            plotlimit = 0.2; figure; hold on;
            plot([x1'; x2'], [y1'; y2'], 'k-');  % faces
            if plotNorms == true
                quiver(face_centers_x, face_centers_y, nx, ny, 0.3, 'r');  % face normals
            end
            scatter(face_centers_x, face_centers_y, 5, 'b', 'filled');   % face centers
            scatter(cell_centers_x, cell_centers_y, 10, 'g', 'filled');  % cell centers
            scatter(ghost_centers(:,1), ghost_centers(:,2), 15, 'm', 'filled'); % ghost cell centers in magenta
            xlim([-plotlimit plotlimit]); ylim([-plotlimit plotlimit]);
            xlabel('x'); ylabel('y');
            title('Face Normals, Face Centers, Cell Centers, and Ghost Cells');
            legend('Faces', 'Normals', 'Face Centers', 'Cell Centers', 'Ghost Cell Centers');
        end
        if plotFill == true
            disp("Moved on to plotting may take awhile heehee")
            rho_vals = Q_numerical(1,:);u_vals = Q_numerical(2,:) ./ rho_vals;   v_vals = Q_numerical(3,:) ./ rho_vals;
            E_vals = Q_numerical(4,:);    p_vals = (gamma - 1) * (E_vals - 0.5 * rho_vals .* (u_vals.^2 + v_vals.^2));
            T_vals = p_vals ./ (rho_vals * 286.69);    velocity_mag = sqrt(u_vals.^2 + v_vals.^2); mach_vals = velocity_mag ./ sqrt(gamma * p_vals ./ rho_vals); 
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
        end
        figure('Name', 'Residual Convergence');semilogy(residual_history, 'LineWidth', 2);
        title('Residual Convergence History');xlabel('Iteration');
        ylabel('Residual (log scale)');legend('L2 Norm Residual');
        if plotCpMach == true
                disp("the mach number and cp on airfoil surface now :D");
                V_inf = sqrt(u_inf^2 + v_inf^2);
                % extract airfoil face centers that we already found 
                airfoil_faces = BC_airfoil_faces;
                x_fc = face_centers(airfoil_faces,1); y_fc = face_centers(airfoil_faces,2);
                n_faces = length(airfoil_faces);
                % normalizing and mapping it to 0 to 1 
                x_le = min(x_fc); x_te = max(x_fc); chord_length = x_te - x_le;x_c = (x_fc - x_le) / chord_length;
                cp_vals = zeros(n_faces,1); mach_vals = zeros(n_faces,1);
                for i = 1:n_faces
                    f = airfoil_faces(i); c = cell_1_of_face(f); % interior cell
                    Q = Q_numerical(:, c); rho = Q(1); u = Q(2)/rho; v = Q(3)/rho; E = Q(4);
                    p = (gamma - 1)*(E - 0.5*rho*(u^2 + v^2));  a = sqrt(gamma * p / rho); V = sqrt(u^2 + v^2);
                    cp_vals(i) = (p - p_inf) / (0.5 * rho_inf * V_inf^2);  mach_vals(i) = V / a;
                end
                is_top = y_fc >= 0; % SPLIT TOP AND BOTTOM BY Y
                is_bot = ~is_top;
                [xtop, itop] = sort(x_c(is_top)); % SORT TOP
                cp_top = cp_vals(is_top); cp_top = cp_top(itop); mach_top = mach_vals(is_top); mach_top = mach_top(itop);
                [xbot, ibot] = sort(x_c(is_bot));  % SORT BOTTOM
                cp_bot = cp_vals(is_bot); cp_bot = cp_bot(ibot); mach_bot = mach_vals(is_bot); mach_bot = mach_bot(ibot);
                % PLOT Cp
                figure; hold on; plot(xtop, cp_top, 'b-', 'LineWidth', 2); plot(xbot, cp_bot, 'r-', 'LineWidth', 2);
                set(gca, 'YDir','reverse'); xlabel('x/c'); ylabel('Cp'); title('Pressure Coefficient on Airfoil'); legend('Top Surface', 'Bottom Surface'); 
                % PLOT Mach
                figure; hold on; plot(xtop, mach_top, 'b-', 'LineWidth', 2); plot(xbot, mach_bot, 'r-', 'LineWidth', 2);
                xlabel('x/c'); ylabel('Mach Number'); title('Mach Number on Airfoil');
                legend('Top Surface', 'Bottom Surface'); grid on;
        end
        if plotContour == true
            disp("smooth contour plots heehee theyre more aesthetic :D");
            % COMPUTE FLOW VARIABLES
            rho_vals = Q_numerical(1,:); u_vals = Q_numerical(2,:) ./ rho_vals; v_vals = Q_numerical(3,:) ./ rho_vals;
            E_vals = Q_numerical(4,:); p_vals = (gamma - 1) * (E_vals - 0.5 * rho_vals .* (u_vals.^2 + v_vals.^2));
            velocity_mag = sqrt(u_vals.^2 + v_vals.^2); mach_vals = velocity_mag ./ sqrt(gamma * p_vals ./ rho_vals);
            % TRIANGULATION BASED ON CELL CENTERS
            tri = delaunay(cell_centers_x, cell_centers_y); tri_centers_x = mean(cell_centers_x(tri), 2); tri_centers_y = mean(cell_centers_y(tri), 2);
            % AIRFOIL POLYGON FROM FACE CENTERS
            airfoil_x = face_centers(BC_airfoil_faces,1); airfoil_y = face_centers(BC_airfoil_faces,2);
            k = convhull(airfoil_x, airfoil_y); airfoil_boundary_x = airfoil_x(k); airfoil_boundary_y = airfoil_y(k);
            % REMOVE TRIANGLES INSIDE AIRFOIL
            inside = inpolygon(tri_centers_x, tri_centers_y, airfoil_boundary_x, airfoil_boundary_y);
            tri = tri(~inside, :);
            % DATA TO PLOT
            contour_vars = {rho_vals, p_vals, velocity_mag, mach_vals};
            contour_titles = {'Density Contour', 'Pressure Contour', 'Velocity Magnitude Contour', 'Mach Number Contour'};
            % LOOP THROUGH SELECTED VARIABLES
            for i = 1:length(contour_vars)
                figure('Name', contour_titles{i});
                trisurf(tri, cell_centers_x, cell_centers_y, contour_vars{i}, 'EdgeColor', 'none', 'FaceColor', 'interp');
                view(2); axis equal tight; colormap jet; colorbar; title(contour_titles{i}); xlabel('x'); ylabel('y');
                xlim([-plotlimit plotlimit]); ylim([-plotlimit plotlimit]); grid off;
            end
        end
        if plotStreamlines == true
            disp("plotting quiver plot because i am the visualization queen :D")
            % COMPUTE FLOW VARIABLES
            rho_vals = Q_numerical(1,:); u_vals = Q_numerical(2,:) ./ rho_vals; v_vals = Q_numerical(3,:) ./ rho_vals; velocity_mag = sqrt(u_vals.^2 + v_vals.^2);
            % GENERATE GRID FOR QUIVERS
            xg = linspace(min(cell_centers_x), max(cell_centers_x), 1000); yg = linspace(min(cell_centers_y), max(cell_centers_y), 3500);
            [X, Y] = meshgrid(xg, yg);  U = griddata(cell_centers_x, cell_centers_y, u_vals, X, Y, 'linear');  V = griddata(cell_centers_x, cell_centers_y, v_vals, X, Y, 'linear');
            % MASK OUT POINTS INSIDE AIRFOIL
            airfoil_x = face_centers(BC_airfoil_faces,1);  airfoil_y = face_centers(BC_airfoil_faces,2);
            k = convhull(airfoil_x, airfoil_y);  mask = inpolygon(X, Y, airfoil_x(k), airfoil_y(k));
            U(mask) = NaN; V(mask) = NaN;
            figure('Name', 'Velocity Magnitude with Quivers'); hold on;
            xlabel('x'); ylabel('y'); title('Velocity Magnitude and Velocity Vectors');
            colormap jet; colorbar; axis equal;
            % PLOT VELOCITY MAGNITUDE
            for i = 1:n_cells
                faces_of_cell = faces_per_cell{i}; xv = [x1(faces_of_cell); x2(faces_of_cell)];
                yv = [y1(faces_of_cell); y2(faces_of_cell)]; k2 = convhull(xv, yv);
                fill(xv(k2), yv(k2), velocity_mag(i), 'EdgeColor', 'none');
            end
            % PLOT QUIVERS SMALL AND BLACK
            quiver(X, Y, U, V, 0.5, 'k'); % 0.5 = scale down arrow length
            xlim([-plotlimit plotlimit]); ylim([-plotlimit plotlimit]);
        end

%% TIME STEPPING FUNCTION
function dt = computeDT(Q_numerical, CFL, Vol, faces_per_cell, nx, ny, lengths, gamma, n_cells, timestep)   
    rho   = Q_numerical(1,:); rho_u = Q_numerical(2,:); rho_v = Q_numerical(3,:); E = Q_numerical(4,:);
    u = rho_u ./ rho; v = rho_v ./ rho; p = (gamma - 1) .* (E - 0.5 .* rho .* (u.^2 + v.^2)); c = sqrt(gamma * p ./ rho);
    dt_local = zeros(n_cells,1);
    for cell = 1:n_cells
        faces = faces_per_cell{cell};
        u_cell = u(cell); v_cell = v(cell); c_cell = c(cell);
        nx_f = nx(faces); ny_f = ny(faces); L_f = lengths(faces); % GET FACE NORMALS AND LENGTHS FOR THIS CELL
        vn = abs(u_cell * nx_f + v_cell * ny_f); % NORMAL VELOCITY MAGNITUDE
        face_sum = sum( (vn + c_cell) .* L_f ); % SUM FACE CONTRIBUTIONS
        dt_local(cell) = CFL * 2 * Vol(cell) / face_sum; % LOCAL TIME STEP
    end
    if timestep == 1
        dt = min(dt_local);
    else
        dt = repmat(dt_local', 4, 1);
    end
end


%% RECONSTRUCTION FUNCTION
function [QL, QR, dQdx, dQdy] = reconstruct(Q, cell_1_of_face, cell_2_of_face, nx, ny, BC_inlet_faces, BC_airfoil_faces, BC_outlet_faces, BC_farfield_faces, ...
                                 mach, gamma, alpha, cell_centers, face_centers, faces_per_cell, order)
    n_faces = length(cell_1_of_face);
    QL = zeros(4, n_faces); QR = zeros(4, n_faces);
    
    rho_inf = 1.4; p_inf = 1; c = sqrt(gamma * p_inf / rho_inf);
    VMAG_inf = c * mach; u_inf = VMAG_inf * cosd(alpha); v_inf = VMAG_inf * sind(alpha);
    E_inf = p_inf / (gamma - 1) + 0.5 * rho_inf * (u_inf^2 + v_inf^2); Q_inf = [rho_inf; rho_inf*u_inf; rho_inf*v_inf; E_inf];
    
    c1 = cell_1_of_face; c2 = cell_2_of_face;
    xf = face_centers(:,1); yf = face_centers(:,2); xc = cell_centers(:,1); yc = cell_centers(:,2);
    
    if order == 1
                %[dQdx, dQdy] = computeGradients(Q, cell_centers, faces_per_cell, cell_1_of_face, cell_2_of_face);

        % FIRST ORDER
        QL = Q(:, c1);
        QR(:, c2 > 0) = Q(:, c2(c2 > 0));
    
        for f = 1:n_faces
            if c2(f) > 0, continue; end
            if ismember(f, BC_inlet_faces)
                QR(:,f) = Q_inf;
            elseif ismember(f, BC_outlet_faces)
                QR(:,f) = Q(:,c1(f));
            elseif ismember(f, BC_farfield_faces)
                QR(:,f) = Q_inf;
                QL(:,f) = Q_inf;
            elseif ismember(f, BC_airfoil_faces)
                qL = Q(:,c1(f));
                rhoL = qL(1); uL = qL(2)/rhoL; vL = qL(3)/rhoL;

                if evalin('base','DIRECT_NS')
                    uR = 0; vR = 0; % no-slip
                else
                    vn = uL * nx(f) + vL * ny(f);
                    uR = uL - 2 * vn * nx(f);  vR = vL - 2 * vn * ny(f);
                end
                QR(:,f) = [rhoL; rhoL*uR; rhoL*vR; qL(4)];
            end
        end

    elseif order == 2
        [dQdx, dQdy] = computeGradients(Q, cell_centers, faces_per_cell, cell_1_of_face, cell_2_of_face);
    
        dx1 = xf - xc(c1); dy1 = yf - yc(c1);
        QL = Q(:,c1) + dQdx(:,c1).*dx1' + dQdy(:,c1).*dy1';
    
        interior = c2 > 0;
        dx2 = xf(interior) - xc(c2(interior)); dy2 = yf(interior) - yc(c2(interior));
        QR(:,interior) = Q(:,c2(interior)) + dQdx(:,c2(interior)).*dx2' + dQdy(:,c2(interior)).*dy2';
    
        for f = 1:n_faces
            if c2(f) > 0, continue; end
            qL = QL(:,f); rhoL = qL(1);  uL = qL(2)/rhoL; vL = qL(3)/rhoL;

            if ismember(f, BC_airfoil_faces)
                if evalin('base','DIRECT_NS')
                    uR = 0; vR = 0; % no-slip
                else
                    vn = uL * nx(f) + vL * ny(f);
                    uR = uL - 2 * vn * nx(f); vR = vL - 2 * vn * ny(f);
                end

                QR(:,f) = [rhoL; rhoL*uR; rhoL*vR; qL(4)];

            elseif ismember(f, BC_inlet_faces) || ismember(f, BC_farfield_faces)
                QR(:,f) = Q_inf;
                if ismember(f, BC_farfield_faces)
                    QL(:,f) = Q_inf;
                end
            elseif ismember(f, BC_outlet_faces)
                QR(:,f) = QL(:,f);
            end
        end
    end
end

function Q_ghost = computeGhostStates(Q, BC_faces, cell_1_of_face, face_centers, ghost_centers, nx, ny)
    n_ghost = length(BC_faces); Q_ghost = zeros(4, n_ghost);
    for k = 1:n_ghost
        f = BC_faces(k);
        c_real = cell_1_of_face(f);
        qR = Q(:,c_real);
    
        rho = qR(1);
        u   = qR(2)/rho;
        v   = qR(3)/rho;
        E   = qR(4);
        % NORMAL DIRECTION
        n = [nx(f), ny(f)];
        vn = u * n(1) + v * n(2);
        % TANGENTIAL COMPONENT STAYS
        u_reflected = u - 2 * vn * n(1);
        v_reflected = v - 2 * vn * n(2);
        % CONSERVED FORM
        Q_ghost(:,k) = [ ...
            rho;
            rho * u_reflected;
            rho * v_reflected;
            E  % no thermal wall assumption yet
        ];
    end
end


function [dQdx, dQdy] = computeGradients(Q, cell_centers, faces_per_cell, cell_1_of_face, cell_2_of_face)
    
    n_cells = size(Q,2); n_var   = size(Q,1); % should be 4
    
    dQdx = zeros(n_var, n_cells); dQdy = zeros(n_var, n_cells);
    
    % ASSUME MAX 10 NEIGHBORS — CHANGE IF NEEDED
    A_max = zeros(10,2); b_max = zeros(10,n_var);
    
    for i = 1:n_cells
        xc = cell_centers(i,1);yc = cell_centers(i,2);
        faces = faces_per_cell{i};n_nb = 0;
    
        for f = faces
            if cell_1_of_face(f) == i
                c2 = cell_2_of_face(f);
            else
                c2 = cell_1_of_face(f);
            end
            if c2 <= 0, continue; end
    
            n_nb = n_nb + 1;
    
            dx = cell_centers(c2,1) - xc; dy = cell_centers(c2,2) - yc;
            dq = (Q(:,c2) - Q(:,i))';
    
            A_max(n_nb, :) = [dx dy]; b_max(n_nb, :) = dq;
        end
    
        % SOLVE IF ENOUGH NEIGHBORS
        if n_nb >= 2
            A_use = A_max(1:n_nb, :);
            b_use = b_max(1:n_nb, :);
    
            grads = lsqminnorm(A_use, b_use); % 2x4 matrix
            dQdx(:,i) = grads(1,:)';
            dQdy(:,i) = grads(2,:)';
        end
    end
end

function Flux = computeFluxes(QL, QR, nx, ny, gamma, dQdx_f, dQdy_f, Fv)

    % COMPUTE INVISCID FLUX (RUSANOV)
    rhoL = QL(1,:); uL = QL(2,:) ./ rhoL; vL = QL(3,:) ./ rhoL; EL = QL(4,:);
    pL = (gamma - 1) .* (EL - 0.5 .* rhoL .* (uL.^2 + vL.^2));
    cL = sqrt(gamma * pL ./ rhoL); vnL = uL .* nx' + vL .* ny';

    rhoR = QR(1,:); uR = QR(2,:) ./ rhoR; vR = QR(3,:) ./ rhoR; ER = QR(4,:);
    pR = (gamma - 1) .* (ER - 0.5 .* rhoR .* (uR.^2 + vR.^2));
    cR = sqrt(gamma * pR ./ rhoR); vnR = uR .* nx' + vR .* ny';

    lambda_max = max([abs(vnL) + cL; abs(vnR) + cR], [], 1);

    FL = [rhoL .* vnL;
          rhoL .* uL .* vnL + pL .* nx';
          rhoL .* vL .* vnL + pL .* ny';
          (EL + pL) .* vnL];

    FR = [rhoR .* vnR;
          rhoR .* uR .* vnR + pR .* nx';
          rhoR .* vR .* vnR + pR .* ny';
          (ER + pR) .* vnR];

    Flux = 0.5 * (FL + FR) - 0.5 * lambda_max .* (QR - QL); % BASE INVISCID FLUX

    % IF NOT DOING NAVIER-STOKES, RETURN
    if ~evalin('base','DIRECT_NS')
        return
    end

    % % SUBTRACT VISCOSITY FROM INTERIOR FACES ONLY
    mask = evalin('base', 'cell_2_of_face') > 0;
    Flux(:, mask) = Flux(:, mask) - Fv(:, mask);
end


function [dQdx_f, dQdy_f] = computeGradientsAtFaces(Q, dQdx, dQdy, ...
    cell_centers, face_centers, cell_1_of_face, cell_2_of_face, ...
    nx, ny, BC_airfoil_faces, wall_type, T_wall)

    [n_var, ~] = size(Q);
    n_faces = length(cell_1_of_face);
    dQdx_f = zeros(n_var, n_faces);
    dQdy_f = zeros(n_var, n_faces);
    
    R = 286.69;
    gamma = 1.4;
    
    for f = 1:n_faces
        i = cell_1_of_face(f);
        j = cell_2_of_face(f);
        n_vec = [nx(f), ny(f)];

        % ====== INTERIOR FACE — PROJECTED BLEND ======
        if j > 0
            ri = cell_centers(i,:); rj = cell_centers(j,:);
            rf = face_centers(f,:);
            li = rf - ri; lj = rf - rj;

            gradQi = [dQdx(:,i), dQdy(:,i)];
            gradQj = [dQdx(:,j), dQdy(:,j)];
            
            dQi = gradQi * li';
            dQj = gradQj * lj';
            
            Qfi = Q(:,i) + dQi;
            Qfj = Q(:,j) + dQj;
            
            Qface = 0.5 * (Qfi + Qfj);
            gradQf = 0.5 * (gradQi + gradQj);
            
            dQdx_f(:,f) = gradQf(:,1);
            dQdy_f(:,f) = gradQf(:,2);

        % ====== WALL FACE — REFLECTED GRADIENT ======
        elseif ismember(f, BC_airfoil_faces)
            qR = Q(:,i); rho = qR(1); u = qR(2)/rho; v = qR(3)/rho; E = qR(4);
            vn = u * n_vec(1) + v * n_vec(2);
            u_g = u - 2 * vn * n_vec(1);
            v_g = v - 2 * vn * n_vec(2);
            qG = [rho; rho*u_g; rho*v_g; E];

            r_real = cell_centers(i,:); r_face = face_centers(f,:);
            proj = dot((r_face - r_real), n_vec);
            r_ghost = r_real + 2 * proj * n_vec;
            dn_vec = r_ghost - r_face; mag_dn = norm(dn_vec);
            if mag_dn < 1e-12, continue; end

            dQ_dn = (qG - qR) / mag_dn;

            % wall energy fix
            if strcmpi(wall_type, 'adiabatic')
                dQ_dn(4) = 0;
            elseif strcmpi(wall_type, 'isothermal')
                V2 = u^2 + v^2;
                E_wall = T_wall * rho * R / (gamma - 1) + 0.5 * rho * V2;
                dQ_dn(4) = (E_wall - qR(4)) / mag_dn;
            end

            dQdx_f(:,f) = dQ_dn * n_vec(1);
            dQdy_f(:,f) = dQ_dn * n_vec(2);
        end
    end

    % ENSURE NO COMPLEX VALUES
    dQdx_f(~isreal(dQdx_f)) = NaN;
    dQdy_f(~isreal(dQdy_f)) = NaN;
end


function Fv = computeViscousFluxes(Q, dQdx_f, dQdy_f, nx, ny, ...
    cell_1_of_face, cell_2_of_face, mu_inf, k_inf, gamma, R)  
    n_faces = size(dQdx_f,2); Fv = zeros(4, n_faces);
    % GET LEFT/RIGHT STATE FOR FACES (LEFT = c1)
    rhoL = Q(1,cell_1_of_face);
    rhoR = Q(1, max(cell_2_of_face,1)); % for boundary faces, repeats left
    
    % GET VELOCITIES (LEFT)
    u = Q(2,cell_1_of_face) ./ rhoL;
    v = Q(3,cell_1_of_face) ./ rhoL;
    
    % DERIVATIVES FOR VELOCITIES
    drho_dx = dQdx_f(1,:);
    drho_dy = dQdy_f(1,:);
    
    du_dx = (dQdx_f(2,:) - u .* drho_dx) ./ rhoL;
    du_dy = (dQdy_f(2,:) - u .* drho_dy) ./ rhoL;
    dv_dx = (dQdx_f(3,:) - v .* drho_dx) ./ rhoL;
    dv_dy = (dQdy_f(3,:) - v .* drho_dy) ./ rhoL;
    
    % DIVERGENCE
    div_v = du_dx + dv_dy;
    
    % STRESS TENSOR (NEWTONIAN, STOKES HYPOTHESIS)
    tau_xx = 2*mu_inf * (du_dx - div_v/3);
    tau_yy = mu_inf * (2*dv_dy - (2/3)*div_v);
    tau_xy = mu_inf * (du_dy + dv_dx);
    
    % ENERGY/EQUATION OF STATE
    E = Q(4,cell_1_of_face);
    V2 = u.^2 + v.^2;
    T = (gamma - 1) * (E - 0.5 * rhoL .* V2) / R;
    
    dT_dx = ((gamma-1)/R) * (dQdx_f(4,:) - u .* dQdx_f(2,:) - v .* dQdx_f(3,:));
    dT_dy = ((gamma-1)/R) * (dQdy_f(4,:) - u .* dQdy_f(2,:) - v .* dQdy_f(3,:));
    qx = -k_inf * dT_dx;
    qy = -k_inf * dT_dy;
    
    % NORMAL PROJECTION
    Fv(1,:) = 0;
    Fv(2,:) = tau_xx .* nx' + tau_xy .* ny';
    Fv(3,:) = tau_xy .* nx' + tau_yy .* ny';
    Fv(4,:) = (u .* tau_xx + v .* tau_xy - qx) .* nx' + ...
              (u .* tau_xy + v .* tau_yy - qy) .* ny';
end


%% RESIDUAL COMPUTATION
function Res = computeResiduals(Flux, cell_1_of_face, cell_2_of_face, lengths, n_cells, Vol)
    Res = zeros(4, n_cells);
    lengths = lengths(:)'; % make row
    % FLUX * LENGTH PER FACE
    Flux_weighted = Flux .* lengths;
    % OUTWARD FLUX FROM CELL 1
    for d = 1:4
        Res(d,:) = Res(d,:) - accumarray(cell_1_of_face, Flux_weighted(d,:)', [n_cells, 1])';
    end
    % INWARD FLUX TO CELL 2 (only valid interior faces)
    mask = cell_2_of_face > 0;
    c2 = cell_2_of_face(mask);
    Flux_in = Flux_weighted(:, mask);
    for d = 1:4
        Res(d,:) = Res(d,:) + accumarray(c2, Flux_in(d,:)', [n_cells, 1])';
    end
    % DIVIDE BY CELL VOLUME
    Res = Res ./ Vol';
end

%% GRID FUNCTIONS
function faces_per_cell = cellsFaces(n_cells, cell_1_of_face, cell_2_of_face)
    faces_per_cell = cell(n_cells, 1);
    for face = 1:length(cell_1_of_face)
        c1 = cell_1_of_face(face); c2 = cell_2_of_face(face);
        if c1 > 0 && c1 <= n_cells
            faces_per_cell{c1}(end+1) = face;
        end
        if c2 > 0 && c2 <= n_cells
            faces_per_cell{c2}(end+1) = face;
        end
    end
end
function Vol = computeCellVolumes(n_cells, faces_per_cell, x1, y1, x2, y2, ~)
    Vol = zeros(n_cells,1);
    for i = 1:n_cells
        faces = faces_per_cell{i};
        xv = [x1(faces); x2(faces)]; yv = [y1(faces); y2(faces)];
        k = convhull(xv, yv); xv = xv(k); yv = yv(k);
        Vol(i) = 0.5 * abs(sum(xv .* yv([2:end 1]) - yv .* xv([2:end 1])));
    end
end
function [airfoil_faces, airfoil_cells] = getAirfoilBoundary(cell_1_of_face, cell_2_of_face, face_centers)
    boundary_faces = find(cell_2_of_face <= 0);
    boundary_x = face_centers(boundary_faces,1);    boundary_y = face_centers(boundary_faces,2);
    is_airfoil = abs(boundary_x) <= 5 & abs(boundary_y) <= 5;
    airfoil_faces = boundary_faces(is_airfoil);    airfoil_cells = cell_1_of_face(airfoil_faces);
end
function [BC_farfield_faces, BC_farfield_cells] = getFarfieldBoundary(cell_1_of_face, cell_2_of_face, BC_inlet_faces, BC_airfoil_faces, BC_outlet_faces)
    is_boundary = (cell_1_of_face <= 0) | (cell_2_of_face <= 0);
    boundary_faces = find(is_boundary);
    used_faces = unique([BC_inlet_faces; BC_airfoil_faces; BC_outlet_faces]);
    BC_farfield_faces = setdiff(boundary_faces, used_faces);
    BC_farfield_cells = zeros(size(BC_farfield_faces));
    for i = 1:length(BC_farfield_faces)
        f = BC_farfield_faces(i);
        if cell_1_of_face(f) > 0
            BC_farfield_cells(i) = cell_1_of_face(f);
        else
            BC_farfield_cells(i) = cell_2_of_face(f);
        end
    end
end
function [inlet_faces, inlet_cells] = getInletBoundary(cell_1_of_face, cell_2_of_face, face_centers)
    boundary_faces = find(cell_2_of_face <= 0);
    boundary_x = face_centers(boundary_faces, 1);
    min_x = min(boundary_x);    tol = 1e-8; 
    inlet_boundary_indices = boundary_faces(abs(boundary_x - min_x) < tol);
    inlet_cells = cell_1_of_face(inlet_boundary_indices);
    inlet_faces = inlet_boundary_indices;
end
function [outlet_faces, outlet_cells] = getOutletBoundary(cell_1_of_face, cell_2_of_face, face_centers)
    boundary_faces = find(cell_2_of_face <= 0);
    boundary_x = face_centers(boundary_faces,1);
    max_x = max(boundary_x);     tol = 1e-6; 
    is_outlet = abs(boundary_x - max_x) <= tol;
    outlet_faces = boundary_faces(is_outlet);
    outlet_cells = cell_1_of_face(outlet_faces);
end
function [n_cells, n_faces, n_points, cell_centers, nx, ny, face_centers, x1, x2, y1, y2, x, y, lengths, cell_1_of_face, cell_2_of_face] ...
          = loadGridData(filename)
    fid = fopen(filename, 'r');
    
    sizes = fscanf(fid, '%d', 3); n_points = sizes(1); n_faces = sizes(2); n_cells = sizes(3); % READ GRID SIZES
    nodes = fscanf(fid, '%f', [2, n_points])'; % NODE COORDINATES
    x = nodes(:, 1); y = nodes(:, 2);
    face_nodes = fscanf(fid, '%d', [2, n_faces])'; % FACE POINTS
    point_1_of_face = face_nodes(:, 1); point_2_of_face = face_nodes(:, 2);
    face_cells = fscanf(fid, '%d', [2, n_faces])'; % FACE-CELL CONNECTIVITY
    cell_1_of_face = face_cells(:, 1); cell_2_of_face = face_cells(:, 2);
    fclose(fid);
    % FACE COORDINATES
    x1 = x(point_1_of_face); y1 = y(point_1_of_face); x2 = x(point_2_of_face); y2 = y(point_2_of_face);
    face_centers = [(x1 + x2)/2, (y1 + y2)/2]; % FACE CENTERS
    
    dx = x2 - x1; dy = y2 - y1; lengths = sqrt(dx.^2 + dy.^2);
    nx = dy ./ lengths; ny = -dx ./ lengths;
    
    % COMPUTE CELL CENTERS (statics way)
    cell_centers = zeros(n_cells, 2);
    face_length_sum = zeros(n_cells, 1);
    
    for f = 1:n_faces
        c1 = cell_1_of_face(f); c2 = cell_2_of_face(f);
        L = lengths(f); fc = face_centers(f,:);
    
        if c1 > 0
            cell_centers(c1,:) = cell_centers(c1,:) + L * fc;
            face_length_sum(c1) = face_length_sum(c1) + L;
        end
        if c2 > 0
            cell_centers(c2,:) = cell_centers(c2,:) + L * fc;
            face_length_sum(c2) = face_length_sum(c2) + L;
        end
    end
    cell_centers = cell_centers ./ face_length_sum;
    % NORMAL DIRECTION TO ALWAYS GO FROM CELL1 TO CELL2
    for f = 1:n_faces
        c1 = cell_1_of_face(f); c2 = cell_2_of_face(f);
        if c1 > 0 && c2 > 0
            cx = cell_centers(c2,1) - cell_centers(c1,1);
            cy = cell_centers(c2,2) - cell_centers(c1,2);
            dot_product = cx * nx(f) + cy * ny(f);
            if dot_product < 0
                nx(f) = -nx(f); ny(f) = -ny(f);
            end
        end
    end
end
function ghost_centers = computeGhostCenters(BC_airfoil_faces, cell_1_of_face, nx, ny, face_centers, cell_centers)
    n_ghost = length(BC_airfoil_faces); ghost_centers = zeros(n_ghost,2);
    for k = 1:n_ghost
        f = BC_airfoil_faces(k);    % face index
        c_real = cell_1_of_face(f); % real cell index
        % VECTOR FROM CELL CENTER TO FACE CENTER
        dx = face_centers(f,1) - cell_centers(c_real,1);
        dy = face_centers(f,2) - cell_centers(c_real,2);
        % REFLECT ACROSS FACE NORMAL
        proj = dx * nx(f) + dy * ny(f);
        xg = cell_centers(c_real,1) + 2 * proj * nx(f);
        yg = cell_centers(c_real,2) + 2 * proj * ny(f);
        ghost_centers(k,:) = [xg yg];
    end
end
function [CD, CL] = computeForceCoefficients(Q_numerical, BC_airfoil_faces, cell_1_of_face, nx, ny, lengths, gamma, rho_inf, V_inf, alpha, chord_length)
    Fx = 0; Fy = 0;
    for f = BC_airfoil_faces'
        c = cell_1_of_face(f); Q = Q_numerical(:,c);
        rho = Q(1); u = Q(2)/rho; v = Q(3)/rho; E = Q(4);
        p = (gamma - 1)*(E - 0.5*rho*(u^2 + v^2)); Fx = Fx + p * nx(f) * lengths(f); Fy = Fy + p * ny(f) * lengths(f);
    end
    q_inf = 0.5 * rho_inf * V_inf^2;
    CD = (Fx * cosd(alpha) + Fy * sind(alpha)) / (q_inf * chord_length); CL = (-Fx * sind(alpha) + Fy * cosd(alpha)) / (q_inf * chord_length);
end
function [dQdx_lim, dQdy_lim] = applyBarthLimiter(Q, dQdx, dQdy, ...
    cell_centers, faces_per_cell, cell_1_of_face, cell_2_of_face)

n_var = size(Q,1);
n_cells = size(Q,2);
dQdx_lim = dQdx; dQdy_lim = dQdy;

for i = 1:n_cells
    xc = cell_centers(i,:); Qi = Q(:,i);
    faces = faces_per_cell{i};
    Qmax = Qi; Qmin = Qi;

    for f = faces
        j = (cell_1_of_face(f) == i) * cell_2_of_face(f) + ...
            (cell_2_of_face(f) == i) * cell_1_of_face(f);

        if j <= 0, continue; end

        xj = cell_centers(j,:);
        dx = xj(1) - xc(1); dy = xj(2) - xc(2);
        Qj = Q(:,j);

        Qmax = max(Qmax, Qj);
        Qmin = min(Qmin, Qj);
    end

    for k = 1:n_var
        dq = dQdx(k,i) * dx + dQdy(k,i) * dy;
        if abs(dq) < 1e-12, continue; end

        r1 = (Qmax(k) - Qi(k)) / dq;
        r2 = (Qmin(k) - Qi(k)) / dq;

        phi = min([1, max(0, min(r1, r2))]);

        dQdx_lim(k,i) = phi * dQdx(k,i);
        dQdy_lim(k,i) = phi * dQdy(k,i);
    end
end
end
