function[Jout] = prec_fun_disk_idx(J,circ_WG_inv,idx_T,idx_R,...
    L_WG,M_WG,N_WG,circ_disk_inv,L_D,M_D,N_D,L,M,N,idxS3,idx_Sq)

% circ_1_inv is the preconditioner for the bus waveguide
% circ_2_inv is the preconditioner for the disk

Jtemp = zeros(3*L*M*N,1);
Jtemp(idxS3) = J;

JtempT = Jtemp(idx_T);
JtempR = Jtemp(idx_R);

% keyboard
JtempT = chan_mvp(circ_WG_inv,JtempT,L_WG,M_WG,N_WG);
% keyboard

JtempR = chan_2_mvp_for_parallel_idx(circ_disk_inv,JtempR,L_D,M_D,N_D,idx_Sq);

Jtemp(idx_T) = JtempT;
%J(idx_B) = JtempB;
Jtemp(idx_R) = JtempR;

Jout = Jtemp(idxS3);
% keyboard