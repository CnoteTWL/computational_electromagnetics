close all; clear all; clc;

cfg = setupParameters();
phi = (0:0.1:359.9) * pi/180; %every angle sweep
% 1 = PEC TMz, 2 = PEC TEz, 3 = dielectric TMz, 4 = dielectric TEz
%for caseID = 1:4
caseID = 1;
    geometry = makeCircleGeometry(cfg.a, cfg.N);
    solution = MoM(geometry, cfg, caseID);
    fields = nearField(solution, geometry, cfg, caseID);
    plotFields(cfg, fields, caseID);
    echo = solveEchoWidth(solution, geometry, cfg, phi);
    Aecho = analyticalEchoWidth(cfg, phi, cfg.N, caseID);
    plotEcho(cfg, echo, Aecho);
    error = sqrt(trapz(cfg.phi, (echo(:)-Aecho(:))).^2) ./ trapz(cfg.phi, Aecho(:));%rmse calculation
    monostatic = interp1(cfg.phi,echo,pi);%monostatic echo width
%end


function cfg = setupParameters()
    cfg.f0 = 10e9;
    u0 = pi*4e-7;
    e0 = 8.85e-12;
    cfg.c0 = 1/sqrt(u0*e0);
    cfg.eta0 = sqrt(u0/e0);
    cfg.a = 0.015;
    cfg.eR = 9;
    cfg.muR = 1;
    cfg.lambda0 = cfg.c0/cfg.f0;
    cfg.lambdaD = cfg.lambda0/sqrt(cfg.eR);
    cfg.k0 = 2*pi/cfg.lambda0;
    cfg.kd = 2*pi/cfg.lambdaD;
    cfg.R = cfg.a+cfg.lambda0;
    cfg.circumference = pi*cfg.lambda0;
    cfg.N = 400;
    cfg.numharmonics = 50;
    cfg.phi = linspace(0,2*pi,361);
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
    geom.yv = yv;
    geom.N = N;
end

function [G, D] = greenStart(x0, y0, geometry, k)
DX = x0(:) - geometry.xmid.';
DY = y0(:) - geometry.ymid.';
R = hypot(DX, DY);
G = (geometry.length.'/(4j)) .* besselh(0, 2, k*R); %green's fn at zero, bessel function is of 2nd kind
D = (k*geometry.length.'/(4j)) .* besselh(1, 2, k*R) .* (geometry.nx.' .* DX + geometry.ny.' .*DY) ./ R; %dirac zero
end

function [G, D] = createGreens(geometry, k)
    [G, D] = greenStart(geometry.xmid, geometry.ymid, geometry, k);
    gamma = exp(0.577);
    d = 1:geometry.N+1:geometry.N^2;
    G(d) = geometry.length/(4j) .* (1 - 1j*(2/pi)*log(gamma*k*geometry.length/(4*exp(1))));
    D(d) = 0;
end

function solution = MoM(geometry, cfg, caseID)
    N = geometry.N;
    I = eye(N); %identity matrix
    Einc = exp(-1j*cfg.k0*geometry.xmid); %incident field of all midpoints
    beta = 1;
    if caseID == 2 || caseID == 4 %checks if TE
        beta = cfg.eR;
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

function echo = solveEchoWidth(solution, geometry, cfg, phi)
cs = cos(phi(:)).';
ss = sin(phi(:)).';
nr = geometry.nx*cs + geometry.ny*ss; %phases needed
%phase = exp(1j*cfg.k0*(geometry.xmid*cs + geometry.ymid*ss)); unnecessary
F = sum(geometry.length .* (solution.q - 1j*nr .* solution.E));
echo = abs(F) .^2 / 4*cfg.k0;
end

function aEcho = analyticalEchoWidth(cfg, phi, N, caseID) %N may be wrong
    n = -N:N;
    x = cfg.k0*cfg.a;
    J = besselj(n,x);
    H = besselh(n, 2, x);
    Js = 0.5*(besselj(n-1,x) - besselj(n+1,x));
    Hs = 0.5*(besselh(n-1,2,x) - besselh(n+1,2,x));
    if caseID == 1 || caseID == 2 %check PEC
        if caseID == 1 || caseID == 3 %check TM
            coeffM = -J ./ H;
        else
            coeffM = -Js ./ Hs;
        end
    else
        xd = cfg.kd*cfg.a; %dielectric case
        Jd = besselj(n, xd);
        Jds = 0.5*(besselj(n-1, xd) - besselj(n+1, xd));
        es = sqrt(cfg.eR);
        if caseID == 2 || caseID == 4 %check TE
            es = 1/es;
        end
        coeffM = -(Js .* Jd - es*J.*Jds) ./ (Hs.*Jd - es*H.*Jds);
    end
    aEcho = (4/cfg.k0) * abs((exp(1j*phi(:)*n) * coeffM.').').^2;
end

function field = nearField(solution, geometry, cfg, caseID)
    x = linspace(-0.06, 0.06, 151);
    [X, Y] = meshgrid(x,x);
    Einc = exp(-1j*cfg.k0*X);
    Escat = zeros(size(X));
    Etot = Escat;
    in = inpolygon(X,Y, geometry.xv, geometry.yv) %returns in or on if in or on edge of polygon
    out = ~in;
    if caseID == 1 || caseID == 3 %check TM
        wavenum = cfg.kd;
    else
        wavenum = cfg.k0;
    end
    [G0, D0] = greenStart(X(out), Y(out), geometry, wavenum);
    Escat(out) = -(G0 * solution.q - D0 * solution.E);
    Etot(out) = Einc(out) + Escat(out);
    field.Et = Etot;
end

function plotEcho(cfg, echo, Aecho)
plot(cfg.phi*180/pi, echo, 'b-', cfg.phi*180/pi, Aecho, 'r--');
xlabel('\phi (deg)'); 
ylabel('Echo width (m)');
legend('MoM', 'Analytical'); 
title('Bistatic Echo Width')
grid on;
end

function plotFields(cfg, fields, caseID)
    imagesc(linspace(-0.06,0.06,size(fields.Et,2)), ...
        linspace(-0.06,0.06,size(fields.Et,1)), abs(fields.Et));
    axis equal tight; 
    cb = colorbar;
    xlabel('x (m)'); ylabel('y (m)');
    title('Field Magnitude');
    if caseID == 1 || caseID == 3
        cb.Label.String = 'Re(E_z) [V/m]';
    else
        cb.Label.String = 'Re(H_z) [A/m]';
    end
end