function[r,idx,dom_x,dom_y,dom_z,xStart,res] = Bragg_geometry(LAM,numPer,width,delW,...
    d,lWG,rWG,len_abs,perturb,res_temp)

dom_x = 2*len_abs + lWG + rWG + numPer*LAM;
dom_y=width+delW;
dom_z = d;

%% Voxelize domain
% Set up preferential direction
h_pref = dom_z;
N_pref = round(h_pref/res_temp);
res = h_pref/N_pref;

[r] = generatedomain_begin_at_x_eq_0(res,dom_x,dom_y/2,dom_z/2);
                   
% plot geometry
xd = r(:,:,:,1);
yd = r(:,:,:,2);
zd = r(:,:,:,3);

xStart = len_abs + lWG;
xEnd = xStart + numPer*LAM;

%% Create geometry for Bragg grating 
% (here we identify the indices of the voxels composing the grating)
idxBraggL = [];
idxBraggR = [];

for i=1:numPer
    x0 = xStart+(i-1)*LAM;
    x1 = xStart+(i-1)*LAM + perturb * res + LAM/2;
    x2 = xStart+i*LAM;
    idxTemp = find((xd>=x0).*(xd<x1));
    idxBraggL = [idxBraggL;idxTemp];
    idxTemp = find((xd>=x1).*(xd<x2).*(abs(yd)<width/2-delW/2));
    idxBraggR = [idxBraggR;idxTemp];
end

idxStart = find((xd<=xStart).*(abs(yd)<=width/2));
idxEnd = find((xd>=xEnd).*(abs(yd)<=width/2));

idx = [idxStart;idxBraggL;idxBraggR;idxEnd];

end