%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% IbIPP - Image-based Initialization and Post-Processing code for 2D %%%
%%%% Topology Optimization, Author: Osezua Ibhadode, October 2020 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ibipp(domain,nelx,volfrac,varargin)
%% Default values of parameters
default_f = 1; default_f_ang = 180; default_pre = 1;
default_pre_load = 0; default_pre_support = 0;
default_rmin = 2; default_filter = 2; default_penal = 3; default_q = 10;
default_beta = 2; default_er = 0.15; default_E0 = 1; default_v = 0.3;
default_tau = 2e-4; default_dist2a = 0; default_theta = 360;
default_opt = 'Density'; default_opt1 = 'SIMP'; default_text = 'none';
default_mtype = 'extrude'; default_exlength = 0.1; default_sym = 'none';

%% Input Parser for required, optional and name-value parameters
pp = inputParser;
addRequired(pp,'domain'); addRequired(pp,'nelx'); addRequired(pp,'volfrac');
addOptional(pp,'forces',default_f); addOptional(pp,'forceAngle',default_f_ang);
addParameter(pp,'optimization',default_opt); addParameter(pp,'densityType',default_opt1);
addParameter(pp,'pressure',default_pre); addParameter(pp,'preserveSupport',default_pre_support);
addParameter(pp,'preserveLoad',default_pre_load); addParameter(pp,'filterRadius',default_rmin);
addParameter(pp,'filter',default_filter); addParameter(pp,'beta',default_beta);
addParameter(pp,'penaltySIMP',default_penal); addParameter(pp,'penaltyRAMP',default_q);
addParameter(pp,'ER',default_er); addParameter(pp,'tau',default_tau);
addParameter(pp,'YoungsModulus',default_E0); addParameter(pp,'PoissonRatio',default_v);
addParameter(pp,'distancetoaxis',default_dist2a); addParameter(pp,'revolutionAngle',default_theta);
addParameter(pp,'modelName',default_text); addParameter(pp,'modelType',default_mtype);
addParameter(pp,'extrudeLength',default_exlength); addParameter(pp,'Symmetry',default_sym);

parse(pp,domain,nelx,volfrac,varargin{:});

%% Optimization parameters
nelx = pp.Results.nelx; volfrac = pp.Results.volfrac;
opt_type = pp.Results.optimization;
IM = pp.Results.densityType; penal = pp.Results.penaltySIMP;
q = pp.Results.penaltyRAMP;
ft = pp.Results.filter; rmin = pp.Results.filterRadius; er = pp.Results.ER;
tau = pp.Results.tau; E0 = pp.Results.YoungsModulus;
v = pp.Results.PoissonRatio;
beta = pp.Results.beta;

%% Load and boundary conditions
Fmag = pp.Results.forces; Fang = pp.Results.forceAngle;
Pmag = pp.Results.pressure;
pre_support = pp.Results.preserveSupport;
pre_load = pp.Results.preserveLoad;

%% Post-processing parameters
tx = pp.Results.modelName; form = pp.Results.modelType;
ht = pp.Results.extrudeLength;
lr = pp.Results.distancetoaxis; theta = pp.Results.revolutionAngle;
symm = pp.Results.Symmetry;

%% Domain Initialization
domain = pp.Results.domain;
[dom,nely] = imageprocessor(domain,nelx);

%% Resolve loads, boundary conditions, non-design domains and preserved regions
[F,fixeddofs,NonD,MusD,volfrac,~] = loadandsupport(nelx,nely,Fmag,Fang,...
    Pmag,dom,pre_support,pre_load,volfrac);

%% Run topology optimization
if strcmpi(opt_type,'Density')
    % SIMP/RAMP density-based
    xPhys = top88mod(nelx,nely,volfrac,penal,rmin,ft,F,fixeddofs,NonD,MusD,...
        IM,q,beta,E0,v);
elseif strcmpi(opt_type,'BESO')
    % Evolutionary
    xPhys = esoLmod(nelx,nely,volfrac,er,rmin,F,fixeddofs,NonD,MusD,E0,v);
elseif strcmpi(opt_type,'levelSet')
    % Level set
    [xPhys,~] = levelset88mod(nelx,nely,volfrac,tau,F,fixeddofs,NonD,MusD,E0,v);
end

%% Post-processing
if ~isempty(find(strcmpi(varargin,tx), 1))
    datatostl(nelx,nely,xPhys,tx,form,ht,lr,theta,symm)
end
end
% This Matlab code was written by Osezua Ibhadode                          %                                     %
% Please sent your comments to: osezua.ibhadode@uwaterloo.ca               %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in a paper which is currently under review                 %  
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but do not guarantee that the code is     %
% free from errors. Furthermore, the author shall not be liable in any     %
% event caused by the use of the program.                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
