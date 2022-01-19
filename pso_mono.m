clc;clear;close all
%% Aircraft Example 
% Monostatic RCS. Freuqnecy 100:150 MHz. 1 degree increment.
% Fine position of scatterers with Particle Swarm Optimization, then find SH coefficients with linear programming.
% Assume real gains

%% Use Precomputed Result
UsePrecomputedLocations = 1;

%% Initial Setting
rng(1)
num_data = 101*121*121; % # of data
K = 16;        % # of scatterers
global f
f = 1e8:1e6:1.5e8;
c = physconst('LightSpeed');
SH_degree = 11;

%% Import rcs Data (Cut out extra data)
addpath('jmontalt-harmonicY-accdfe4')
rcs = importdata('rcsData/1degree_100MHz200MHz_subangle.txt');
f_data = rcs.data(:,1);
f_data = reshape(f_data, [101, 121, 121]);
f_data = f_data(1:size(f, 2), :, :);

phi_data = rcs.data(:,2);
phi_data = reshape(phi_data, [101, 121, 121]);
phi_data = phi_data(1:size(f, 2), :, :);

theta_data = rcs.data(:,3);
theta_data = reshape(theta_data, [101, 121, 121]);
theta_data = theta_data(1:size(f, 2), :, :);

rcs_Re = rcs.data(:,6)*2*sqrt(pi); % 2*sqrt(pi) is from a CST offset 
rcs_Re = reshape(rcs_Re, [101, 121, 121]);
rcs_Re = rcs_Re(1:size(f, 2), :, :);

rcs_Im = rcs.data(:,7)*2*sqrt(pi); % 2*sqrt(pi) is from a CST offset 
rcs_Im = reshape(rcs_Im, [101, 121, 121]);
rcs_Im = rcs_Im(1:size(f, 2), :, :);

figure
for n=1:size(f, 2)
    plotImage(30:150, 210:330, squeeze(pow2db(rcs_Re(n,:,:).^2+rcs_Im(n,:,:).^2)), append('RCS from CST at ',num2str((n-1)+100),'MHz (dB)'))
    pause(0.03)
end

%% Construct Real Spherical Harmonics
L = (SH_degree+1)^2;
SH_matrix = [];
[PHI, THETA] = ndgrid(210:330, 30:150);
for l=0:SH_degree
    for m=-l:l
        SH_matrix = cat(3, SH_matrix, harmonicY(l, m, deg2rad(THETA), deg2rad(PHI), 'type', 'real'));
%         % Visualize spherical harmonics
%         figure
%         [x, y, z] = sph2cart(deg2rad(PHI), pi/2-deg2rad(THETA), abs(SH_matrix));
%         surf(x(:,:,end), y(:,:,end), z(:,:,end), SH_matrix(:,:,end));
%         axis equal
%         xlabel('X')
%         ylabel('Y')
%         zlabel('Z')
    end
end
global H
H = reshape(SH_matrix, [121*121, L]);

%% Construct a list of direction unit vectors with theta_data and phi_data
rho = 1;
r_xy = rho .* sind(reshape(THETA, [121*121,1]));
x_dir  = r_xy  .* cosd(reshape(PHI, [121*121,1]));
y_dir  = r_xy  .* sind(reshape(PHI, [121*121,1]));
z_dir  = rho .* cosd(reshape(THETA, [121*121,1]));
global direction
direction = [x_dir y_dir z_dir];

%% Use PSO to find locations
rcs_Re = reshape(permute(rcs_Re, [2,3,1]), [121*121, size(f, 2)]);
rcs_Im = reshape(permute(rcs_Im, [2,3,1]), [121*121, size(f, 2)]);
global res_Re
global res_Im
res_Re = rcs_Re;
res_Im = rcs_Im;
location_all =[];
S_bar_all = [];
r_bar = [reshape(res_Re, [size(res_Re,1)* size(res_Re,2), 1]); reshape(res_Im, [size(res_Im,1)* size(res_Im,2), 1])];
disp(['Initial residual norm: ', num2str(norm(r_bar))])
if UsePrecomputedLocations==0
    tic
    for k=1:K
        location = (particleswarm(@evalLocation2, 3, [-10,-10,-10],[10,10,10]))';
        toc
        warning('off','all')
        location_all = [location_all location];
        lead = 2*direction*location/c;
        S_temp = exp(1i*2*pi*repelem(f', 121*121).*repmat(lead, size(f, 2), 1));
        S = repmat(H, size(f, 2), 1).*S_temp;
        S_bar_all = [S_bar_all, [real(S); imag(S)]];       
%         res = r_bar - S_bar_all*(S_bar_all\r_bar);
        res = r_bar - S_bar_all*((S_bar_all'*S_bar_all + diag(repelem(0.001,L*k)))\(S_bar_all'*r_bar));
        res_Re = reshape(res(1:size(res, 1)/2), [121*121, size(f, 2)]);
        res_Im = reshape(res(size(res, 1)/2+1:end), [121*121,size(f, 2)]);
        disp([num2str(k),'th ', 'scatterer found. ', 'New residual norm: ', num2str(norm(res))])
    end
    toc
    disp('All scatterers found')
%     alpha = S_bar_all\r_bar;
    alpha = (S_bar_all'*S_bar_all + diag(repelem(0.001,L*K)))\(S_bar_all'*r_bar);
    rcs_hat = S_bar_all*alpha;
    rcs_hat = rcs_hat(1:size(res, 1)/2) + 1i*rcs_hat(size(res, 1)/2+1:end);
    % save('location_all.mat', 'location_all')
    % save('alpha.mat', 'alpha')
    % save('rcs_hat.mat', 'rcs_hat')
else
    addpath('results')
    load('location_all.mat');
    load('alpha.mat');
    load('rcs_hat.mat')
end
rcs_hat_real_mat = reshape(real(rcs_hat), [121,121,size(f, 2)]);
rcs_hat_imag_mat = reshape(imag(rcs_hat), [121,121,size(f, 2)]);
rcs_Re_mat = reshape(rcs_Re, [121, 121, size(f, 2)]);
rcs_Im_mat = reshape(rcs_Im, [121, 121, size(f, 2)]);

%% Plot Comparisons Animation
for n=1:51
    subplot(2,2,1)
    plotImage(30:150, 210:330, rcs_Re_mat(:,:,n), append('RCS from CST at ',num2str((n-1)+100),'MHz'))
    caxis([-40 40])
    subplot(2,2,2)
    plotImage(30:150, 210:330, rcs_Im_mat(:,:,n), append('RCS from CST at ',num2str((n-1)+100),'MHz'))
    caxis([-40 40])
    subplot(2,2,3)
    plotImage(30:150, 210:330, squeeze(rcs_hat_real_mat(:,:,n)), append('Modeled RCS at ',num2str((n-1)+100),'MHz'))
    caxis([-40 40])
    subplot(2,2,4)
    plotImage(30:150, 210:330, squeeze(rcs_hat_imag_mat(:,:,n)), append('Modeled RCS at ',num2str((n-1)+100),'MHz'))
    caxis([-40 40])
    pause(0.01)
end

%% Plot Comparisons varying f
figure
hold on;
set(gca,'FontSize',10)
plot(f'/1e6,squeeze(rcs_Re_mat(1,1,:)),'b','LineWidth',2)
plot(f'/1e6,squeeze(rcs_hat_real_mat(1,1,:)),'b--','LineWidth',2)
plot(f'/1e6,squeeze(rcs_Re_mat(61,71,:)),'r','LineWidth',2)
plot(f'/1e6,squeeze(rcs_hat_real_mat(61,71,:)),'r--','LineWidth',2)
xlabel('Signal frequency (MHz)') 
ylabel('RCS value') 
plot(f'/1e6,squeeze(rcs_Im_mat(1,1,:)),'g','LineWidth',2)
plot(f'/1e6,squeeze(rcs_hat_imag_mat(1,1,:)),'g--','LineWidth',2)
plot(f'/1e6,squeeze(rcs_Im_mat(61,71,:)),'m','LineWidth',2)
plot(f'/1e6,squeeze(rcs_hat_imag_mat(61,71,:)),'m--','LineWidth',2)
legend('CST (\theta = 30^o, \phi = 210^o, Re)','Model (\theta = 30^o, \phi = 210^o, Re)',...
    'CST (\theta = 90^o, \phi = 280^o, Re)','Model (\theta = 90^o, \phi = 280^o, Re)',...
    'CST (\theta = 30^o, \phi = 210^o, Im)','Model (\theta = 30^o, \phi = 210^o, Im)',...
    'CST (\theta = 90^o, \phi = 280^o, Im)','Model (\theta = 90^o, \phi = 280^o, Im)')

%% Plot RCS Comparisons from CST at 100MHz
figure
subplot(1,2,1)
plotImage(30:150, 210:330, rcs_Re_mat(:,:,1), append('RCS from CST at 100MHz (Real)'))
caxis([-40 40])
subplot(1,2,2)
plotImage(30:150, 210:330, rcs_Im_mat(:,:,1), append('RCS from CST at 100MHz (Imaginary)'))
caxis([-40 40])
figure
subplot(1,2,1)
plotImage(30:150, 210:330, squeeze(rcs_hat_real_mat(:,:,1)), append('Modeled RCS at 100MHz (Real)'))
caxis([-40 40])
subplot(1,2,2)
plotImage(30:150, 210:330, squeeze(rcs_hat_imag_mat(:,:,1)), append('Modeled RCS at 100MHz (Imaginary)'))
caxis([-40 40])

%% Plot Comparisons varying theta
figure
hold on;
set(gca,'FontSize',10)
plot(30:150,squeeze(rcs_Re_mat(:,61,1)),'b','LineWidth',2)
plot(30:150,squeeze(rcs_hat_real_mat(:,61,1)),'b--','LineWidth',2)
xlabel('\theta') 
ylabel('RCS value') 
plot(30:150,squeeze(rcs_Im_mat(:,61,1)),'g','LineWidth',2)
plot(30:150,squeeze(rcs_hat_imag_mat(:,61,1)),'g--','LineWidth',2)
legend('CST (\phi = 270^o, Re)','Model (\phi = 270^o, Re)',...
    'CST (\phi = 270^o, Im)','Model (\phi = 270^o, Im)')

%% Plot Point Scatterer Model
SH_matrix = [];
[PHI, THETA] = ndgrid(0:1:360, 0:1:180);
for l=0:SH_degree
    for m=-l:l
        SH_matrix = cat(3, SH_matrix, harmonicY(l, m, deg2rad(THETA), deg2rad(PHI), 'type', 'real'));
    end
end
X = sind(THETA).*cosd(PHI);
Y = sind(THETA).*sind(PHI);
Z = cosd(THETA);
% reflectionGain_max = [];
% reflectionGain_min = [];
figure
for k=1:K
    X_k = X/1.2 + location_all(1, k);
    Y_k = Y/1.2 + location_all(2, k);
    Z_k = Z/1.2 + location_all(3, k);
%     reflectionGain_max = [reflectionGain_max max(abs(sum(SH_matrix .* X_hat_mat_, 3)),[],'all')];
%     reflectionGain_min = [reflectionGain_min min(abs(sum(SH_matrix .* X_hat_mat_, 3)),[],'all')];
    surf(X_k, Y_k, Z_k, sum(SH_matrix .* reshape(alpha(1+L*(k-1):L*k), [1,1,L]), 3))
    hold on
end
axis equal
colorbar
shading interp
caxis('auto')
% disp([ 'Max abs of reflection gains: ',num2str(max(reflectionGain_max))]);
% disp([ 'Max abs of reflection gains: ',num2str(min(reflectionGain_min))]);

%% Plot Residual and Computation Time
r = [6432.6351; 5419.4704; 4621.9973; 4060.8951; 3612.2708; 3306.5446; 3072.5149; 2825.9404; 2693.2167; 2565.7404; 2485.6842; 2371.0728; 2269.1932; 2183.8941; 2110.1319; 2051.6988; 1983.0872];
t = [0;1338.829591; 2930.124026; 3929.019936; 5122.144224; 7741.889500;  8943.516396; 10179.553625; 11483.886351; 13094.897263; 14162.679922; 15631.833317; 16701.522530; 17889.725334; 19612.064301; 20769.547664; 22513.469068];
figure
set(gca,'FontSize',15)
yyaxis left
plot(0:16,r,'--o','LineWidth',2)
xlabel('Number of scatterers') 
ylabel('Residual') 
yyaxis right
plot(0:16, t/60/60,'--o','LineWidth',2) 
ylabel('Computation Time (Hr)') 

%% Plot 1 scatterer residual varying of scatterer location in x,y,z (has to be run before the algorithm)
figure
cost_all=[];
cost_all2=[];
cost_all3=[];
for x = -10:0.1:10
cost = evalLocation([x,0,0]); cost_all=[cost_all cost];
end
for y = -10:0.1:10
cost = evalLocation([0,y,0]); cost_all2=[cost_all2 cost];
end
for z = -10:0.1:10
cost = evalLocation([0,0,z]); cost_all3=[cost_all3 cost];
end
hold on
plot(-10:0.1:10, cost_all, 'LineWidth',2)
plot(-10:0.1:10, cost_all2, 'LineWidth',2)
plot(-10:0.1:10, cost_all3, 'LineWidth',2)
set(gca,'FontSize',10)
legend('x','y','z')
xlabel('Scatterer position') 
ylabel('Cost') 