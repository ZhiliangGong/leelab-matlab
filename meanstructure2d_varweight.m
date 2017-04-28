function [meanstruct, bestframe, trj, vel] = meanstructure2d_varweight(trj, trjbf, index, mass, tolerance, vel)
%% meanstructure2d
% calc average structure by iterative superimpose in the XY space (Z-axis ignored)
%
%% Syntax
%# crd = meanstructure2d(trj);
%# crd = meanstructure2d(trj);
%# [crd, trj] = meanstructure2d(trj);
%# [crd, trj] = meanstructure2d(trj, index_atom);
%# [crd, trj] = meanstructure2d(trj, index_atom, mass);
%# [crd, trj] = meanstructure2d(trj, index_atom, mass, tolerance);
%# [crd, trj] = meanstructure2d(trj, index_atom, mass, tolerance, vel);
%# [crd, trj, vel] = meanstructure2d(trj, [], [], [], vel);
%
%% Description
% This routine calculates the average structure from given
% trajectory. The algorithm superimposes the trajectories to a
% plausible average structure, then updates the average structrue.
% In the superimpose step, only XY-space is considered (Z-direction
% is ignored, this style should be convenient for membrane proteins).
% This process is repeated until some convergence is achieved in
% the RMSD between the average structures.
% This routine may be useful for a preprocess
% for the subsequent structure-analysis routines, such as Principal
% Component Analysis. 
%
%% Example
%# trj = readnetcdf('ak.nc');
%# [crd, trj] = meanstructure(trj);
%
%% See also
% superimpose
%

%% initialization
ref = trj(1, :);
weight = mass;
rmsd = realmax;

if ~exist('index', 'var')
  index = [];
end
  
if ~exist('mass', 'var')
  mass = [];
end

if ~exist('tolerance', 'var') || isempty(tolerance)
  tolerance = 10^(-6);
end

if ~exist('vel', 'var')
  vel = [];
end

%% iterative superimpose
ref = decenter2d(ref, index, mass);
trj = decenter2d(trj, index, mass);


while rmsd > tolerance
  weight_old = weight;
  ref_old = ref;
  [~, trj, vel] = superimpose2d(ref_old, trj, index, weight_old, vel, true);
  ref = mean(trj);
  weight = (mass')./(std(trj(:,1:3:end-2)).^2 + std(trj(:,2:3:end-1)).^2 + std(trj(:,3:3:end)).^2);
  weight = weight';
  trj = decenter2d(trj, index, weight);
  rmsd = superimpose2d(ref_old, ref, index, weight_old, [], true);
  fprintf('rmsd from the previous mean structure: %f A\n', rmsd);
end

if numel(vel) ~= 0
  vel = decenter2d(vel, index, weight_old);
end

meanstruct.crd = ref;
meanstruct.crdstd = mass./weight;
[~, trj, vel] = superimpose2d(ref, trj, index, weight, vel);

%% find representative frame that best fits mean structure

trj_3d = zeros(size(trjbf));
[ref_3d, cow] = decenter(ref, index, weight);
trj_3d(:,1:3:end-2) = trjbf(:,1:3:end-2) - cow(1);
trj_3d(:,2:3:end-1) = trjbf(:,2:3:end-1) - cow(2);
trj_3d(:,3:3:end) = trjbf(:,3:3:end) - cow(3);
rmsd_3d = zeros(size(trj_3d,1), 1);
for i = 1:size(trj,1)
    rmsd_3d(i) = sum(superimpose(trj_3d(i,:), ref_3d, index, weight, [], true));
end

[bestframe.rmsd, imin] = min(rmsd_3d);
fprintf('rmsd best frame, %i, from mean structure: %f A\n', imin, bestframe.rmsd);
bestframe.struct = trj_3d(imin,:);
bestframe.framenum = imin-1;

