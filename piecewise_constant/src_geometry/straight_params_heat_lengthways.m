function[Mr,Mc] = ...
    straight_params_heat_lengthways(free_space_lambda,n_ext,material,len_abs,...
    r,dom_x,idx,e_r_left,e_r_rate,lWG,rWG)

% free_space_lambda : Free-space wavelength of incident field
% n_ext             : refractive index of exterior medium, e.g, Air, SiO2.
% material          : String telling us the material
% len_abs           : absorber length
% r                 : coordinates of voxel centers
% dom_x             : length of domain in x-direction (longest dimension)

mu = 4*pi*1e-7;
co = 299792458;             % speed of light in vacuum
freq = co./free_space_lambda * n_ext; 
eo = 1/co^2/mu;
omega = 2*pi*freq;

[L,M,N,~] = size(r);

if strcmp(material,'Si')==1
    e_r_free = Lorentz(free_space_lambda);
    n_r = sqrt(e_r_free)/n_ext;    % scale by exterior medium
    e_r = n_r^2;
end

%% Set up absorbers
xd = r(:,:,:,1);
z_abs_R = dom_x-len_abs;      % location of beginning of right absorber
z_abs_L  = len_abs;  % location of beginning of left absorber

% fun_absorb = @(z) (z>=z_abs_R).*(abs(z-z_abs_R)./len_abs).^2 +...
%     (z<=z_abs_L).*(abs(z-z_abs_L)./len_abs).^2;

abs_profile=2;
factors = [1 1/2 1/3 1/4 1/5 0.403653];
factor = factors(abs_profile+1);
% keyboard

R0=1e-10;
e_r_right = e_r;
sigma0_left = -log(R0)/len_abs * 1/factor*1/4 * ...
                sqrt(eo*e_r_left/mu) * 2; % factor 2 disappears for matched impendances

sigma0_right = -log(R0)/len_abs * 1/factor*1/4 * ...
                sqrt(eo*e_r_right/mu) * 2;
            
%    keyboard         
            
%% Varying the temperature in the long direction of the grating
x_Bragg_start = len_abs + lWG;
x_Bragg_end = dom_x - len_abs - rWG;
fun_perm = @(x) (x<=x_Bragg_start).*e_r_left + ...
                (x>x_Bragg_start).*(x<x_Bragg_end).*(e_r_left + ...
                (x-x_Bragg_start).*e_r_rate) + ...
                (x>=x_Bragg_end).*e_r_right;
            
 
fun_absorb = @(z) (z>=z_abs_R).*sigma0_right.*(abs(z-z_abs_R)./len_abs).^2 +...
    (z<=z_abs_L).*sigma0_left.*(abs(z-z_abs_L)./len_abs).^2;            
    

%%

sigma_e = zeros(L,M,N);
sigma_e(idx) = fun_absorb(xd(idx));

epsilon_r = ones(L,M,N);
epsilon_r(idx) = fun_perm(xd(idx)); 

epsilon_r = epsilon_r - 1j*sigma_e/(eo*omega);

Mr = epsilon_r;
Mc = (epsilon_r - 1.0);

