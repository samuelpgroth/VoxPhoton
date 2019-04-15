%% Si_in_SiO2_parameters.m
% Establish frequencies and refractive indices of Si and surrounding SiO2
% For the frequency sweep, we use wavelengths located at Legendre points in
% the 1500nm-1600nm range

n_ext = 1.444;
s_plus = 1600e-9;
s_minus = 1500e-9;
s_mid = (s_plus+s_minus)/2;

% Legengre points
[wG,xG]=gauss_1d(nG);

s = s_mid + (s_plus-s_minus)/2 * xG;
free_wavelength = (s_minus:1e-9:s_plus);
numLambda = length(free_wavelength);
co = 299792458;             % speed of light in vacuum
free_freq = co./free_wavelength * 1.444; 
freq_plus = co/s_plus * 1.444;
freq_minus = co/s_minus * 1.444;
freq_s = co./s .* 1.444;