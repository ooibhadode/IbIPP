# IbIPP
 Image-based Initialization and Post-Processing code for Topology Optimization
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%% IbIPP - Image-based Initialization and Post-Processing code for 2D %%%
 %%%% Topology Optimization. Author: Osezua Ibhadode October, 2020.%%%%%%%%%
 %    IbIPP is an open-source code for initializing and post-processing
 %   free-form 2D topology optimization problems.
 %
 %   IBIPP(FILE, NELX, VOLFFRAC) Performs topology optimization on the image
 %   provided in FILE using the number of elements in the x-direction, NELX,
 %   the required volume fraction, VOLFRAC and other default parameters.
 %
 %   IBIPP(FILE, NELX, VOLFRAC, FMAG, FANG) Performs topology optimization
 %   using FILE, NELX, VOLFRAC, vectors of force magnitudes, FMAG, and
 %   force angles, FANG, and other default parameters
 %
 %   IBIPP(...,'PropertyName',VALUE,'PropertyName',VALUE,...) performs
 %   topology optimization using the following property values:
 %.............................................................................
 %    Name-value        Description
 %..............................................................................
 %   pressure        .  Pressure loads (multiple loads must be in a vector)
 %   optimization    .  Topology optimization approach  
 %   densityType     .  Type of density-based method
 %   preserveSupport .  Preserve elements within the support region in the image
 %   preserveLoad    .  Preserve elements within the load region in the image
 %   filterRadius    .  Specify filter radius for density-based and BESO approaches
 %   filter          .  Specify filter type for density-based approaches
 %   beta            .  Regularization parameter for Heaviside projection
 %   penaltySIMP     .  Penalty value for SIMP density-based approach
 %   penaltyRAMP     .  Penalization factor for RAMP density-based approach
 %   ER              .  Evolution ratio for BESO
 %   tau             .  Regularization parameter for level set method
 %                   .  (reaction diffusion)
 %   YoungsModulus   .  Youngs Modulus material property
 %   PoissonRatio    .  Poisson Ratio material property
 %   modelName       .  Name given to model obtained by extrusion or
 %                   .  revolution of the optimized topology. Must end in .stl
 %   modelType       .  Modeling by extrusion or revolution
 %   extrudeLength   .  Length for extrusion. It is a factor multiplied by
 %                   .  the minimum of nelx and nely to obtain the length of
 %                   .  the 3rd dimension
 %   symmetry        .  Specifies the position of the line of symmetry for
 %                   .  symmetry-based modeling
 %   distancetoaxis  .  Distance between the optimized topology and its
 %                   .  center of revolution to the left for revolution-based modeling
 %   revolutionAngle .  The angle of revolution for revolution-based modeling
 %
 % %   Table: Color codes (RGB) for creating an input image file
 %   -----------------------------------------------------------------------
 %   Features                Color code rules            Recommended values                                       
 %   -----------------------------------------------------------------------
 %   Designable domain       R∩G∩B<200                      0     0     0
 %   Non-designable domain	R∩G∩B>200                     255   255   255
 %   Point load              200≤R*≤255,G=0,B=0            200    0     0
 %   Pressure                200≤R≤230,100≤G*≤150,B=0     200   100    0
 %   Preserved region        R=0,200≤G≤255,B=0              0    200    0
 %   Fixed region            R=0,G=0,200≤B≤255              0     0    200
 %   Region fixed along x	100≤R≤150,G=0,200≤B≤255      100    0    200
 %   Region fixed along y	R=0,200≤G≤255,200≤B≤255       0    200   200        
 %   * for multiple forces or pressure, R should be in multiples of 5
 %   starting from 200 for forces (200, 205, 210,...etc) and G should be
 %   in multiples of 5 starting from 100 for pressures (100, 105, 110,..etc)
 %
 %   Example 1:
 %     % Optimize a half-MBB using default values of optional inputs
 %    ibipp('ex1.png',150,0.4)
 %
 %   Example 2:
 %     % Optimize a Hammerhead pier using the BESO approach
 %     ibipp('ex2.png',200,0.5,[2 1 1 2],[180 180 180 180],'optimization',...
 %   'BESO','filterRadius',3)
 %
 %   Example 3:
 %     % Optimize a 2-point loading mechanical part while preserving load
 %     and support elements
 %     ibipp('ex3.png',300,0.45,[1 1],[0 0],'optimization','levelset',...
 %   'tau',5e-5,'preserveLoad',1,'preserveSupport',1)
 %
 %   Example 4:
 %     % Optimize a half-spanner design and generate an extruded full
 %     spanner STL model
 %    ibipp('half_spanner.png',250,0.5,'pressure',[1,1],'preserveload',2,...
 %   'preservesupport',1,'modelname','spanner.stl','symmetry','left',...
 %   'modeltype','extrude','extrudelength',0.2)
 %
 %   IbIPP utilizes other open-source codes such as top88.m by
 %   Andreassen et.al, esoL.m by Xia et. al., levelset88.m by Otomori et.al,
 %   revolve2D.m by Treeby and Cox, and stlwrite.m by Sven.
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
