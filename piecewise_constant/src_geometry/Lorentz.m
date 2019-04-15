function[ep] = Lorentz(lambda)

epsilon = 7.9874;
ep_lorentz = 3.6880;
omega0 = 3.9328*1e15;
c = 299792458;

ep = epsilon + ep_lorentz*omega0^2/(omega0^2 - (2*pi*c/lambda)^2);