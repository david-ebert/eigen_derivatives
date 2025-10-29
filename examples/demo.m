% demonstrates the functions of this package
%
% PLEASE START VIA EXAMPLE FILE
% See also EXAMPLE_NONDEGENERATE, EXAMPLE_CROSS_FRISWELL, EXAMPLE_DEFLECT,
% EXAMPLE_CROSS_CROSS, EXAMPLE_CROSS_DEFLECT, EXAMPLE_PAIR_CROSS

addpath ..

%%

x_ = -1:.01:1;
y_ = -1:.01:1;

n_x = length(x_);
n_y = length(y_);

t_ = linspace(-1,1,101);
n_t = length(t_);

%% build parameterized Matrix

K_fris =@(x,y) K_0 ...
    + x*dxK + x.^2/2*dx2K + x.^3/6*dx3K + x.^4/24*dx4K...
    + y*dyK + y.^2/2*dy2K + y.^3/6*dy3K + y.^4/24*dy4K ;

%% EVP in reference point

K0 = K_fris(x0,y0);
[u0,l0] = eig(K0,"vector");
l0 = l0(1:n_ev);
u0 = u0(:,1:n_ev);

%% directional derivatives of straight path

% direction
r = 0; % in [0,1];
% r = rand;

angle_x = r*2*pi + pi/4;
dx = cos(angle_x); dy = sin(angle_x);

dK_{1} = dx         *dxK  + dy         *dyK ...
       + dx*x0      *dx2K + dy*y0      *dy2K ...
       + dx*x0^2/2  *dx3K + dy*y0^2/2  *dy3K ...
       + dx*x0^3/6  *dx4K + dy*y0^3/6  *dy4K;
dK_{2} = dx^2       *dx2K + dy^2       *dy2K ...
       + dx^2*x0    *dx3K + dy^2*y0    *dy3K ...
       + dx^2*x0^2/2*dx4K + dy^2*y0^2/2*dy4K;
dK_{3} = dx^3       *dx3K + dy^3       *dy3K ...
       + dx^3*x0    *dx4K + dy^3*y0    *dy4K;
dK_{4} = dx^4       *dx4K + dy^4       *dy4K;

%% derivatives of eigenpair

% derivatives (with respect to eigenspace)
tic
[dL_,du_] = eig_der(l0(1),u0,K0,dK_);
disp("Derivatives (with respect to eigenspace) in " + num2str(toc) + " seconds.")

if ~degenerate
    dl_ = dL_;
else
    % polarization
    tic
    [P0,k_pol] = pol(dL_);
    disp("Initial polarization in " + num2str(toc) + " seconds.")
    disp('Polarization:')
    disp(k_pol)
    % derivatives of polarization
    tic
    [dP_,dl_] = pol_der(dL_,P0,k_pol,u0,du_);
    disp("Derivatives of polarization in " + num2str(toc) + " seconds.")
    % polarized derivatives of eigenfunctions
    [du_pol,u0_pol] = du_pol(u0,du_,P0,dP_);
end

%% print derivatives of eigenvalues

disp('Eigenvalue derivatives:')
tab = table(l0,dl_{:});
tab.Properties.VariableNames = {'l0' 'dl' 'd2l' 'd3l' 'd4l'};
disp(tab);

%% Taylor approximations

%#ok<*SAGROW>
for i_t = 1:n_t
    t = t_(i_t);

    l_Taylor{1}(:,i_t) = l0                  + t*     dl_{1};
    l_Taylor{2}(:,i_t) = l_Taylor{1}(:,i_t)  + t^2/2* dl_{2};
    l_Taylor{3}(:,i_t) = l_Taylor{2}(:,i_t)  + t^3/6* dl_{3};
    l_Taylor{4}(:,i_t) = l_Taylor{3}(:,i_t)  + t^4/24*dl_{4};

    if ~degenerate
        u_Taylor{1,i_t} = u0              + t*     du_{1};
        u_Taylor{2,i_t} = u_Taylor{1,i_t} + t^2/2* du_{2};
        u_Taylor{3,i_t} = u_Taylor{2,i_t} + t^3/6* du_{3};
    else
        u_Taylor{1,i_t} = u0_pol          + t*     du_pol{1};
        u_Taylor{2,i_t} = u_Taylor{1,i_t} + t^2/2* du_pol{2};
        u_Taylor{3,i_t} = u_Taylor{2,i_t} + t^3/6* du_pol{3};
    end
end

%% evaluation for surface plot

l_surf = zeros(n_inner,n_x,n_y);
for i_x = 1:n_x
    for i_y = 1:n_y
        l_surf(:,i_x,i_y) = eig(K_fris(x_(i_x),y_(i_y)),"vector");
    end
end

%% plot

colors = colororder;

figure(1); clf
subplot(1,2,1); 
hold on
plot3(x0,y0,l0,'k*:')

x_traj = x0+dx*t_;
y_traj = y0+dy*t_;

for i_ev = 1:n_inner
    surf(x_,y_,squeeze(l_surf(i_ev,:,:))')
end
for i_ev = 1:n_ev
    % plot3(x_traj,y_traj,l_Taylor{1}(i_ev,:),'.:','Color',colors(i_ev,:),'DisplayName','Taylor1')
    % plot3(x_traj,y_traj,l_Taylor{2}(i_ev,:),'.-.','Color',colors(i_ev,:),'DisplayName','Taylor2')
    % plot3(x_traj,y_traj,l_Taylor{3}(i_ev,:),'.--','Color',colors(i_ev,:),'DisplayName','Taylor3')
    plot3(x_traj,y_traj,l_Taylor{4}(i_ev,:),'.-','Color',colors(i_ev,:),'DisplayName','Taylor4')
end
view(3)
alpha(.9)
shading interp
xlabel('x'); ylabel('y'); zlabel('\lambda')
title('Surfaces')
xlim([-1,1]);ylim([-1,1]);
zlim([1,3]); clim([1,3])
colorbar

subplot(1,2,2);
hold on
plot(0,l0,'k*:')
for i_ev = 1:n_ev
    % plot(t_,l_Taylor{1}(i_ev,:),'.:','Color',colors(i_ev,:),'DisplayName','Taylor1')
    % plot(t_,l_Taylor{2}(i_ev,:),'.-.','Color',colors(i_ev,:),'DisplayName','Taylor2')
    % plot(t_,l_Taylor{3}(i_ev,:),'.--','Color',colors(i_ev,:),'DisplayName','Taylor3')
    plot(t_,l_Taylor{4}(i_ev,:),'.-','Color',colors(i_ev,:),'DisplayName','Taylor4')
end
xlim([-1,1])
ylim([1,3])
xlabel('t'); ylabel('\lambda')
title('Trajectory')