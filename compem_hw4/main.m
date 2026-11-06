close all; clear all; clc;

cfg = setupParameters();
phi = (0:0.1:359.9) * pi/180; %every angle sweep
% 1 = PEC TMz, 2 = PEC TEz, 3 = dielectric TMz, 4 = dielectric TEz
for caseID = 1:4
    
end


function cfg = setupParameters()
    cfg.f0 = 10e9;
    cfg.c0 = 1/sqrt(pi*4e-7*8.85e-12);
    cfg.a = 0.015;
    cfg.epsilonR = 9;
    cfg.muR = 1;
    cfg.lambda0 = cfg.c0/cfg.f0;
    cfg.lambdaD = cfg.lambda0/sqrt(epsilonR);
    cfg.k0 = 2*pi/cfg.lambda0;
    cfg.R = cfg.a+cfg.lambda0;
    cfg.circumference = pi*lambda0;
    cfg.numplot = 400;
    cfg.numharmonics = 50;
    % cfg.hmax = cfg.lambda0/30;
    % cfg.phi = linspace(0,2*pi,361);
    % cfg.Rc = cfg.a+cfg.lambda0/4; % Far field integration contour
    % cfg.Nc = 1440;
end

function makeCircleGeom(a,N)
    phv = 2*pi*(0:N) .' /N; %angle vector 0-2pi
    xv = a*cos(phv);
    yv = a*sin(phv);
    dx = diff(xv); %finds differences of adjacent components
    dy = diff(yv);
    geom.length = hypot(dx, dy); %hypotenuse of sides dx, dy
    geom.xmid = (xv(1:end-1) + xv(2:end))/2; %returns midpoint vector of x vector
    geom.ymid = (yv(1:end-1) + yv(2:end))/2;
    geom.nx = dy./geom.length; %calculates x component of unit normal vector
    geom.ny = -dx./geom.length;
    geom.xv = xv; %geom. here is unnecessary but more compact
    geom.yv = yh;
    geom.N = N;
end

function MoM(geometry, cfg, caseID)
    N = geometry.N;
    I = eye(N); %identity matrix
    Einc = exp
end

%G->J
%D->M