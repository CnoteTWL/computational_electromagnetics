close all; clear all; clc;

cfg = setupParameters();
phi = (0:0.1:359.9) * pi/180; %every angle sweep
% 1 = PEC TMz, 2 = PEC TEz, 3 = dielectric TMz, 4 = dielectric TEz
for caseID = 1:4
    geometry = makeCircleGeometry(cfg.a, cfg.N);
    solution = MoM(geometry, cfg, caseID);
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
    cfg.N = 400;
    cfg.numharmonics = 50;
    % cfg.hmax = cfg.lambda0/30;
    % cfg.phi = linspace(0,2*pi,361);
    % cfg.Rc = cfg.a+cfg.lambda0/4; % Far field integration contour
    % cfg.Nc = 1440;
end

function geom = makeCircleGeometry(a,N)
    phv = 2*pi*(0:N) .' /N; %angle vector evenly spaced 0-2pi
    xv = a*cos(phv); %x and y coordinates of the vertices
    yv = a*sin(phv);
    dx = diff(xv); %finds differences of adjacent vertices
    dy = diff(yv);
    geom.length = hypot(dx, dy); %hypotenuse of sides dx, dy
    geom.xmid = (xv(1:end-1) + xv(2:end))/2; %returns midpoint vector of x vector
    geom.ymid = (yv(1:end-1) + yv(2:end))/2;
    geom.nx = dy./geom.length; %calculates x component of outward unit normal vector
    geom.ny = -dx./geom.length;
    geom.xv = xv; %geom. here is unnecessary but more compact
    geom.yv = yh;
    geom.N = N;
end

function [G, D] = greenStart(x0, y0, geometry, k)
DX = x0(:) - geometry.xm.';
DY = y0(:) - geometry.ym.';
R = hypot(DX, DY);
G = (geometry.L.'/(4j)) .* besselh(0, 2, k*R); %green's fn at zero, bessel function is of 2nd kind
D = (k*geometry.L.'/(4j)) .* besselh(1, 2, k*R) .* (geometry.nx.' .* DX + geometry.ny.' .*DY) ./ R; %dirac zero
end

function [G, D] = createGreens(geometry, k)
    [G, D] = greenStart(geometry.xmid, geometry.ymid, geometry, k);
    gamma = exp(0.577);
    d = 1:geometry.N+1:geometry.N^2;
    G(d) = geometry.L/(4j) .* (1 - 1j*(2/pi)*log(gamma*k*geometry.L/(4*exp(1))));
    D(d) = 0;
end

function solution = MoM(geometry, cfg, caseID)
    N = geometry.N;
    I = eye(N); %identity matrix
    Einc = exp(-1j*cfg.k0*geometry.xmid); %incident field of all midpoints
    beta = 1;
    if caseID == 2 || caseID == 4 %checks if TE
        beta = cfg.epsr;
    end
    [G, D] = createGreens(geometry, cfg.k0);
    if caseID == 1 || caseID == 2 %checks if PEC
        if caseID == 1 || caseID == 3 %checks if TM
            q = G \ Einc;
            E = zeros(N,1);
        else
            E = (0.5*I - D) \ Einc;
            q = zeros(N,1);
        end
    end
    A = [0.5*I - D, G;
        0.5*I + D, -beta.*G];
    x = A \ [Einc; zeros(N,1)];
    E = x(1:N);
    q = x(N+1:end);
    %wrap outputs
    solution.E = E;
    solution.q = q;
    solution.beta = beta;
    %define surface currents
    if caseID == 1 || caseID == 3
        solution.J_z = q / (1j*cfg.k0*cfg.eta0);
    end
    if caseID == 2 || caseID == 4
        solution.J_t = -E;
    end
end

