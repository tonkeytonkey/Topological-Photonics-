function [omega,fem] = SquareRodTMWGpec(kx)

flclear fem

% COMSOL version
clear vrsn
vrsn.name = 'COMSOL 3.4';
vrsn.ext = '';
vrsn.major = 0;
vrsn.build = 248;
vrsn.rcs = '$Name:  $';
vrsn.date = '$Date: 2007/10/10 16:07:51 $';
fem.version = vrsn;

% Geometry
g3=rect2('1','15.0','base','center','pos',{'0','0'},'rot','0');
g6=circ2('0.175','base','center','pos',{'0','-7.0'},'rot','0');
%g6=ellip2(0.24,0.12,'base','center','pos',[0,-4],'rot','60');
garr=geomarrayr(g6,0,1,1,15);



% Analyzed geometry
clear s
s.objs={g3,g6,garr{:}};


fem.draw=struct('s',s);
fem.geom=geomcsg(fem);

% Initialize mesh
fem.mesh=meshinit(fem, ...
                  'hauto',3);

% Constants
fem.const = {'kx',kx, ...
  'eps1','15', ...
  'mu','0.84', ...
  'kr','0.41*i'};

% (Default values are not included)

% Application mode 1
clear appl
appl.mode.class = 'InPlaneWaves';
appl.module = 'RF';
appl.sshape = 2;
appl.assignsuffix = '_rfwe';
clear prop
prop.analysis='eigen';
prop.eigtype='lambda';
appl.prop = prop;
clear bnd
bnd.type = {'E0','cont','periodic'};
bnd.index = {1,0,0};
bnd.kper = {{'kx';0},{0;0},{'kx';0}};
bnd.pertype = {'floque','sym','floque'};
bnd.ind = [3,1,1,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, ...
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2];
appl.bnd = bnd;
clear equ
equ.epsilonr = {1,'eps1'};
equ.mur = {1,{'mu','kr',0;'-kr','mu',0;0,0,'mu'}};
equ.ind = [1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2];
appl.equ = equ;
fem.appl{1} = appl;
fem.border = 1;

% Multiphysics
fem=multiphysics(fem);

% Extend mesh
fem.xmesh=meshextend(fem);

% Solve problem
fem.sol=femeig(fem, ...
               'complexfun','on', ...
               'conjugate','on', ...
               'solcomp',{'Ez'}, ...
               'outcomp',{'Ez'}, ...
               'neigs',350, ...
               'shift',-8e6i);   %-i*0.47*2*pi*3e8


omega = -imag(fem.sol.lambda)/2/pi/3e8;
end