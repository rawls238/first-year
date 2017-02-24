load('simu_calvo_opt.mat')
[~,~,Vy_opt,Vx_opt] = variance_decomposition(gx,hx,nETASHOCK);
Vx_sum = sum(Vx_opt);
Vy_sum = sum(Vy_opt);
Vxr_yT = [Vx_opt(1,nyT)/Vx_sum(nyT); Vx_opt(2,nyT)/Vx_sum(nyT)]
Vyr_ht = [Vy_opt(1,nh)/Vy_sum(nh); Vy_opt(2,nh)/Vy_sum(nh)]
Vyr_ct = [Vy_opt(1,ncT)/Vy_sum(ncT); Vy_opt(2,ncT)/Vy_sum(ncT)]
Vxr_pt = [Vx_opt(1,np)/Vx_sum(np); Vx_opt(2,np)/Vx_sum(np)]
Vxy_tby = [Vy_opt(1,ntby)/Vy_sum(ntby); Vy_opt(2,ntby)/Vy_sum(ntby)]

load('simu_calvo_peg.mat')
[~,~,Vy_peg,Vx_peg] = variance_decomposition(gx,hx,nETASHOCK);
Vx_sum = sum(Vx_peg);
Vy_sum = sum(Vy_peg);
Vxr_yT = [Vx_peg(1,nyT)/Vx_sum(nyT); Vx_peg(2,nyT)/Vx_sum(nyT)]
Vyr_ht = [Vy_peg(1,nh)/Vy_sum(nh); Vy_peg(2,nh)/Vy_sum(nh)]
Vyr_ct = [Vy_peg(1,ncT)/Vy_sum(ncT); Vy_peg(2,ncT)/Vy_sum(ncT)]
Vxr_pt = [Vx_peg(1,np)/Vx_sum(np); Vx_peg(2,np)/Vx_sum(np)]
Vxy_tby = [Vy_peg(1,ntby)/Vy_sum(ntby); Vy_peg(2,ntby)/Vy_sum(ntby)]
    

