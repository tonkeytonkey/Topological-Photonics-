
close all; 
clear all; 
clc; 

f = @(x,y) exp(cos(sqrt(x.^2 + y.^2))); % a function of two variables
d = -2*pi:0.1:2*pi;                     % domain for both x,y
[X,Y] = meshgrid(d,d);                  % create a grid of points
Z = f(X,Y);                             % evaluate f at every point on grid

[nrows,ncols] = size(X);                % first obtain the size of X
Z1 = f(X(:),Y(:));                      % convert X,Y to column vectors and evaluate f
Z1 = reshape(Z1,nrows,ncols);           % reshape

f4 = figure;                            % create a new figure
p9 = surf(X,Y,Z);                       % plot the surface of the function
shading interp;                         % interpolate between the points
material dull;                          % alter the reflectance
%camlight(90,0);                         % add some light - see doc camlight
%alpha(0.8);                             % make slightly transparent
box on;  