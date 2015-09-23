function [omega,fem] = SquareRodTM(bands,kx,ky,rod,epsilon,mu,khr)

%{
bands=1;
kx=0;
ky=0;
rod=0.3;
epsilon=2;
mu=1;
khr=0; 
  %}
     
flclear fem;

% COMSOL version
clear vrsn
vrsn.name = 'COMSOL 3.5';
vrsn.ext = 'a';
vrsn.major = 0;
vrsn.build = 603;
vrsn.rcs = '$Name:  $';
vrsn.date = '$Date: 2008/12/03 17:02:19 $';
fem.version = vrsn;

% Constants
fem.const = {
  'kx',kx, ...
  'ky',ky};

% Geometry
g1=square2('1','base','corner','pos',{'0','0'},'rot','0');
%g2=square2('0.5','base','corner','pos',{'0','0'},'rot','0');
g2=circ2(rod,'base','center','pos',{'0.5','0.5'},'rot','0');
%g2=rect2('0.12','0.12','base','center','pos',[0.5,0.5],'rot','0','const',fem.const);
%g2=square2(0.12,'base','center','pos',{'0.5','0.5'},'rot','0');
%g2=ellip2(0.24,0.12,'base','center','pos',[0.5,0.5],'rot','45','const',fem.const);

% Analyzed geometry
clear s;
s.objs={g1,g2};
s.name={'SQ1','R1'};
s.tags={'g1','g2'};

fem.draw=struct('s',s);
fem.geom=geomcsg(fem);

% Initialize mesh
fem.mesh=meshinit(fem, ...
                  'hauto',3);
              
% Application mode 1
clear appl
appl.mode.class = 'InPlaneWaves';
appl.module = 'RF';
appl.assignsuffix = '_rfwe';
clear prop
prop.analysis='eigen';
prop.field='TE'; % comsol TE, our TM
appl.prop = prop;
clear bnd
bnd.pertype = {'floque','sym','floque'};
bnd.type = {'periodic','cont','periodic'};
bnd.index = {2,0,1};
bnd.kper = {{'kx';'ky'},{0;0},{'kx';'ky'}};
bnd.ind = [3,1,1,3,2,2,2,2];
appl.bnd = bnd;
clear equ
equ.epsilonr = {1,{1;1;'eps1'}};
equ.mur = {1,{'mur1','i*kr',0;'-i*kr','mur2',0;0,0,1}};
%equ.mur = {1,{'mur1','kr',0;'kr','mur2',0;0,0,1}};
equ.ind = [1,2];
appl.equ = equ;
fem.appl{1} = appl;
fem.sdim = {'x','y'};
fem.frame = {'ref'};
fem.border = 1;
clear units;
units.basesystem = 'SI';
fem.units = units;

% Scalar expressions


fem.expr = {
  'eps1',epsilon, ...
  'mur1',mu, ...
  'mur2','mur1', ...
  'kr',khr};
%descr.expr= {'mur','Relative permeability tensor element','gamma','Gyromagnetic ratio','wm','Larmor frequency at saturation limit','H0','Applied magnetic bias field','Ms','Saturation magnetization','kr','Relative permeability tensor element','w0','Larmor frequency'};
%fem.descr = descr;
% % Scalar expressions
% fem.expr = {'gamma','0*1.759e11[C/kg]', ...
%   'H0','omega_rfwe/(gamma*mu0_rfwe+1e4[m/C])', ...
%   'w0','mu0_rfwe*gamma*H0', ...
%   'Ms','0.3[Wb/m^2]/mu0_rfwe', ...
%   'wm','mu0_rfwe*gamma*Ms', ...
%   'mur','1+w0*wm/(w0^2-omega_rfwe^2)', ...
%   'kr','omega_rfwe*wm/(w0^2-omega_rfwe^2)'};
% descr.expr= {'mur','Relative permeability tensor element','gamma','Gyromagnetic ratio','wm','Larmor frequency at saturation limit','H0','Applied magnetic bias field','Ms','Saturation magnetization','kr','Relative permeability tensor element','w0','Larmor frequency'};
% fem.descr = descr;

% ODE Settings
clear ode
clear units;
units.basesystem = 'SI';
ode.units = units;
fem.ode=ode;
% Multiphysics
fem=multiphysics(fem);

% Extend mesh
fem.xmesh=meshextend(fem);

Nsolv = bands;
% Solve problem
fem.sol=femeig(fem, ...
               'complexfun','on', ...
               'conjugate','on', ...
               'solcomp',{'Ez'}, ...
               'outcomp',{'Ez'}, ...
               'neigs',Nsolv ...%               'eigref','4.5e8'
               );
           
omega = -imag(fem.sol.lambda)/2/pi/3e8;
omega = omega(1:Nsolv);
end