function[ R ] = chimapplot(R, rot)
%Function to create plot of chi^2 values of R structure output of
%Ref_analysis of protein orientations

theta_rotcount = 0;
theta = zeros((rot(3)-rot(1))/rot(2) + 1, 1);
phi = zeros((rot(6)-rot(4))/rot(5) + 1, 1);
chi = zeros(((rot(6)-rot(4))/rot(5)) + 1, (rot(3)-rot(1))/rot(2) + 1);

for theta_rot = rot(1):rot(2):rot(3)
    theta_rotcount = theta_rotcount + 1;
    phi_rotcount = 0;
    theta(theta_rotcount) = theta_rot;
    for phi_rot = rot(4):rot(5):rot(6)
        phi_rotcount = phi_rotcount + 1;
        chi(phi_rotcount, theta_rotcount) = R.fitdata.(sprintf('t%03dp%03d', theta_rot, phi_rot)).chi;
        if theta_rot == rot(1)
            phi(phi_rotcount) = phi_rot; %only compute phi array once
        end
    end
end

[Chimin, ~] = min(chi(:));

pmap = [0,0,0.562500000000000;0,0,0.598958313465118;0,0,0.635416686534882;...
    0,0,0.671875000000000;0,0,0.708333313465118;0,0,0.744791686534882;...
    0,0,0.781250000000000;0,0,0.817708313465118;0,0,0.854166686534882;...
    0,0,0.890625000000000;0,0,0.927083313465118;0,0,0.963541686534882;...
    0,0,1;0,0.200000002980232,1;0,0.400000005960465,1;0,0.600000023841858,1;...
    0,0.800000011920929,1;0,1,1;0.142857149243355,1,0.857142865657806;...
    0.285714298486710,1,0.714285731315613;0.428571432828903,1,0.571428596973419;...
    0.571428596973419,1,0.428571432828903;0.714285731315613,1,0.285714298486710;...
    0.857142865657806,1,0.142857149243355;1,1,0;1,0.946708440780640,0;...
    1,0.893416941165924,0;1,0.840125381946564,0;1,0.786833882331848,0;...
    1,0.733542323112488,0;1,0.680250763893127,0;1,0.626959264278412,0;...
    1,0.573667705059052,0;1,0.520376205444336,0;1,0.467084646224976,0;...
    1,0.413793116807938,0;1,0.331034481525421,0;1,0.248275876045227,0;...
    1,0.165517240762711,0;1,0.0827586203813553,0;1,0,0;0.978260874748230,0,0;...
    0.956521749496460,0,0;0.934782624244690,0,0;0.913043498992920,0,0;...
    0.891304373741150,0,0;0.869565188884735,0,0;0.847826063632965,0,0;...
    0.826086938381195,0,0;0.804347813129425,0,0;0.782608687877655,0,0;...
    0.760869562625885,0,0;0.739130437374115,0,0;0.717391312122345,0,0;...
    0.695652186870575,0,0;0.673913061618805,0,0;0.652173936367035,0,0;...
    0.630434811115265,0,0;0.608695626258850,0,0;0.586956501007080,0,0;...
    0.565217375755310,0,0;0.543478250503540,0,0;0.521739125251770,0,0;...
    0.500000000000000,0,0];

dof = length(R.fitdata.refnorm(:,1))-6-1;

R.chimap.theta = theta;
R.chimap.phi = phi;
R.chimap.chiarray = chi;
R.chimap.chinormarray = chi/Chimin;
chisym = sym(dof*R.chimap.chiarray);
R.chimap.likelihoodarray = vpa(exp(-(1/2)*chisym));
R.chimap.likelihoodarray = R.chimap.likelihoodarray/sum(R.chimap.likelihoodarray(:));
R.chimap.likelihoodarray = double(R.chimap.likelihoodarray);
chipvalue = fcdf(R.chimap.chinormarray, dof, dof);
R.chimap.pvaluearray = ones(size(chipvalue)) - chipvalue;
[Iphimin, Ithetamin] = find(chi == Chimin);
if length(Ithetamin) == 1
    R.chimap.min = R.fitdata.(sprintf('t%03dp%03d', theta(Ithetamin), phi(Iphimin)));
    R.chimap.min.theta = theta(Ithetamin);
    R.chimap.min.phi = phi(Iphimin);
else
    for i = 1:length(Ithetamin)
        R.chimap.(['min', num2str(i)]) = R.fitdata.(sprintf('t%03dp%03d', theta(Ithetamin(i)), phi(Iphimin(i))));
        R.chimap.(['min', num2str(i)]).theta = theta(Ithetamin(i));
        R.chimap.(['min', num2str(i)]).phi = phi(Iphimin(i));
    end
end

CorrFact = length(R.chimap.theta)*length(R.chimap.phi);

figure(1)
imagesc(R.chimap.theta,R.chimap.phi,R.chimap.chinormarray);
set(gca,'FontSize', 15, 'FontName', 'Helvetica'), set(gca,'YDir','normal');
ylabel('\phi (deg)'); xlabel('\theta (deg)'); 
colorbar;
colormap(pmap); 
caxis([1 3]); 
colorbar('YTick',[finv(0.95,dof,dof), finv(0.99,dof,dof), finv(0.999,dof,dof), finv(0.9999,dof,dof), finv(0.99999,dof,dof)],...
    'YTickLabel',{'P=0.05','.01','.001','1e-4','<1e-5'});
title('\chi^2 Map by F test P Value', 'FontSize', 20)

figure(2)
imagesc(R.chimap.theta,R.chimap.phi,R.chimap.chinormarray);
set(gca,'FontSize', 15, 'FontName', 'Helvetica'), set(gca,'YDir','normal');
ylabel('\phi (deg)'); xlabel('\theta (deg)'); 
colorbar;
colormap(pmap); 
caxis([1 7])
colorbar('YTick',[finv((1-0.05)^(1/(CorrFact)),dof,dof), finv((1-0.01)^(1/(CorrFact)),dof,dof), finv((1-0.001)^(1/(CorrFact)),dof,dof), finv((1-0.0001)^(1/(CorrFact)),dof,dof), finv((1-0.00001)^(1/(CorrFact)),dof,dof)],...
    'YTickLabel',{'P=0.05','.01','.001','1e-4','<1e-5'});
title('\chi^2 Map by F test P Value Corrected', 'FontSize', 20)

figure(3)
imagesc(R.chimap.theta,R.chimap.phi,(R.chimap.chiarray - Chimin*ones(size(R.chimap.chiarray))));
set(gca,'FontSize', 15, 'FontName', 'Helvetica'), set(gca,'YDir','normal');
ylabel('\phi (deg)'); xlabel('\theta (deg)'); 
colorbar;
colormap(pmap);
caxis([0, chi2inv(0.999999,2)]);
colorbar('YTick',[chi2inv(0.95,2), chi2inv(0.99,2), chi2inv(0.999,2), chi2inv(0.9999,2), chi2inv(0.99999,2)],...
    'YTickLabel',{'P=0.05','.01','.001','1e-4','<1e-5'});
title('\chi^2 Map by \chi^2 P Value', 'FontSize', 20)

figure(4)
surf(R.chimap.theta,R.chimap.phi,R.chimap.likelihoodarray);
set(gca,'FontSize', 15, 'FontName', 'Helvetica'), set(gca,'YDir','normal');
ylabel('\phi (deg)'); xlabel('\theta (deg)'); 
colorbar;
colormap(flipud(pmap));
title('Likelihood Surface Plot', 'FontSize', 20);

end