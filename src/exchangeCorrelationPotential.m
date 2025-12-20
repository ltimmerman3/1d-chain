function [S] = exchangeCorrelationPotential(S)
rho = S.rho;
rho(rho < S.xc_rhotol) = S.xc_rhotol;
if S.isgradient
    drho = S.grad * rho;
    sigma = drho.*drho;
    sigma(sigma < S.xc_rhotol) = S.xc_rhotol;
end

% iexch 
switch S.ixc(1)
    case 1
        [ex,vx] = slater(rho);
        v2x = zeros(size(rho));
    case 2
        [ex,vx,v2x] = pbex(rho,sigma,S.xc_option(1));
    otherwise
        ex = zeros(size(rho));
        vx = zeros(size(rho));
        v2x = zeros(size(rho));
end

% icorr
switch S.ixc(2)
    case 1
        [ec,vc] = pz(rho);
        v2c = zeros(size(rho));
    case 2
        [ec,vc] = pw(rho);
        v2c = zeros(size(rho));
    case 3
        [ec,vc,v2c] = pbec(rho,sigma,S.xc_option(2));
    otherwise
        ec = zeros(size(rho));
        vc = zeros(size(rho));
        v2c = zeros(size(rho));
end  

if S.usefock > 1 % hybrid
    if S.xc == 427
        [exsr,v1xsr,v2xsr] = pbexsr(rho,sigma,S.hyb_range_pbe);
        ex = ex - S.hyb_mixing * exsr ./ rho;
        vx = vx - S.hyb_mixing * v1xsr;
        v2x = v2x - S.hyb_mixing * v2xsr;
    else
        ex = ex * (1-S.hyb_mixing);
        vx = vx * (1-S.hyb_mixing);
        v2x = v2x * (1-S.hyb_mixing);
    end
end

exc = ex + ec;
vxc = vx + vc;
v2xc = v2x + v2c;

if S.isgradient % xc involves gradient of rho
    vxc = vxc - S.grad * (v2xc.*drho);
end

S.e_xc = exc;
S.Vxc = vxc;
S.dvxcdgrho = v2xc;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ex,vx] = slater(rho)
% slater exchange
% @param rho  = total electron density

% parameter 
C2 = 0.73855876638202; % 3/4*(3/pi)^(1/3)
C3 = 0.9847450218427;  % (3/pi)^(1/3)

% computation
ex = - C2 * rho.^(1./3.);
vx = - C3 * rho.^(1./3.);
end

function [ec,vc] = pw(rho)
% pw correlation
% @param rho  = total electron density
% J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)

% parameters
A = 0.031091 ;
alpha1 = 0.21370 ;
beta1 = 7.5957 ;
beta2 = 3.5876 ;
beta3 = 1.6382 ;
beta4 = 0.49294 ;

% computation
rs = (0.75./(pi*rho)).^(1./3.);
rsm12 = rs.^(-0.5);
rs12 = rs.^0.5;
rs32 = rs.^1.5;
rs2 = rs.^2;

om = 2*A*(beta1*rs12 + beta2*rs + beta3*rs32 + beta4*rs2);
dom = A*(beta1*rsm12+ 2*beta2 + 3*beta3*rs12 + 2*2*beta4*rs);
olog = log(1 + 1./om);
t = -2*A*(1+alpha1*rs);
ec = t.*olog;
vc = ec - (rs/3.).*(-2*A*alpha1*olog - (t.*dom)./(om.*om+om) ) ;
end

function [ec,vc] = pz(rho)
% pz correlation
% @param rho  = total electron density
% J.P. Perdew and A. Zunger, PRB 23, 5048 (1981).

% parameters
A = 0.0311;
B = -0.048;
C = 0.002;
D = -0.0116;
gamma1 = -0.1423;
beta1 = 1.0529;
beta2 = 0.3334;

% compuatation
ec = zeros(size(rho,1),1);
vc = zeros(size(rho,1),1);
rs = (0.75./(pi*rho)).^(1.0/3.0) ;
islt1 = (rs < 1.0);
lnrs = log(rs(islt1));
sqrtrs = sqrt(rs(~islt1));
ec(islt1) = A * lnrs + B + C * rs(islt1) .* lnrs + D * rs(islt1);
ox = 1.0 + beta1*sqrtrs + beta2*rs(~islt1);
ec(~islt1) = gamma1 ./ ox;
vc(islt1) = lnrs.*(A + (2.0/3.0)*C*rs(islt1)) + (B-(1.0/3.0)*A) + (1.0/3.0)*(2.0*D-C)* rs(islt1);
vc(~islt1) = ec(~islt1) .* (1 + (7.0/6.0)*beta1*sqrtrs + (4.0/3.0)*beta2*rs(~islt1)) ./ ox;
end

function [ex,v1x,v2x] = pbex(rho,sigma,iflag)
% pbe exchange:
% @param rho  = total electron density
% @param grho = |\nabla rho|^2
% @param iflag options
% iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
% iflag=2  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
% iflag=3  RPBE: B. Hammer, et al., Phys. Rev. B 59, 7413 (1999)
% iflag=4  Zhang-Yang Revised PBE: Y. Zhang and W. Yang., Phys. Rev. Lett. 80, 890 (1998)
assert(iflag == 1 || iflag == 2 || iflag == 3 || iflag == 4);

% parameters 
mu_ = [0.2195149727645171 10.0/81.0 0.2195149727645171  0.2195149727645171];
mu = mu_(iflag);
kappa_ = [0.804 0.804 0.804 1.245];
kappa = kappa_(iflag);
threefourth_divpi = 3.0/4.0/pi;
sixpi2_1_3 = (6.0 * pi^2)^(1.0/3.0);
sixpi2m1_3 = 1.0/sixpi2_1_3;
mu_divkappa = mu/kappa;

% computation
rho_updn = rho/2.0;
rho_updnm1_3 = rho_updn.^(-1.0/3.0);
rhomot = rho_updnm1_3;
ex_lsd = -threefourth_divpi * sixpi2_1_3 * (rhomot .* rhomot .* rho_updn);
rho_inv = rhomot .* rhomot .* rhomot;
coeffss = (1.0/4.0) * sixpi2m1_3 * sixpi2m1_3 * (rho_inv .* rho_inv .* rhomot .* rhomot);
ss = (sigma/4.0) .* coeffss;

if iflag == 1 || iflag == 2 || iflag == 4
    divss = 1.0./(1.0 + mu_divkappa * ss);
    dfxdss = mu * (divss .* divss);
elseif iflag == 3
    divss = exp(-mu_divkappa * ss);
    dfxdss = mu * divss;
end

fx = 1.0 + kappa * (1.0 - divss);
dssdn = (-8.0/3.0) * (ss .* rho_inv);
dfxdn = dfxdss .* dssdn;
dssdg = 2.0 * coeffss;
dfxdg = dfxdss .* dssdg;

ex = ex_lsd .* fx;
v1x = ex_lsd .* ((4.0/3.0) * fx + rho_updn .* dfxdn);
v2x = 0.5 * ex_lsd .* rho_updn .* dfxdg;
end

function [ec,v1c,v2c] = pbec(rho,sigma,iflag)
% pbe correlation 
% @param rho  = total electron density
% @param sigma = |\nabla rho|^2
% iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
% iflag=2  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
% iflag=3  RPBE: B. Hammer, et al., Phys. Rev. B 59, 7413 (1999)
assert(iflag == 1 || iflag == 2 || iflag == 3);

% parameter 
beta_ = [0.066725 0.046 0.066725];
beta = beta_(iflag);
rsfac = 0.6203504908994000; % (0.75/pi)^(1/3)
sq_rsfac = sqrt(rsfac);
sq_rsfac_inv = 1.0/sq_rsfac;
third = 1.0/3.0;
twom1_3 = 2.0^(-third);
ec0_aa = 0.031091; 
ec0_a1 = 0.21370;  
ec0_b1 = 7.5957;
ec0_b2 = 3.5876;   
ec0_b3 = 1.6382;   
ec0_b4 = 0.49294;  
gamma = (1.0 - log(2.0)) /pi^2;
gamma_inv = 1/gamma;
coeff_tt = 1.0/(4.0 * 4.0 / pi * (3.0 * pi^2)^third);

% computation
rho_updn = rho/2.0;
rho_updnm1_3 = rho_updn.^(-third);
rhom1_3 = twom1_3 * rho_updnm1_3;

rhotot_inv = rhom1_3 .* rhom1_3 .* rhom1_3;
rhotmo6 = sqrt(rhom1_3);
rhoto6 = rho .* rhom1_3 .* rhom1_3 .* rhotmo6;

rs = rsfac * rhom1_3;
sqr_rs = sq_rsfac * rhotmo6;
rsm1_2 = sq_rsfac_inv * rhoto6;

%        Formulas A6-A8 of PW92LSD
ec0_q0 = -2.0 * ec0_aa * (1.0 + ec0_a1 * rs);
ec0_q1 = 2.0 * ec0_aa *(ec0_b1 * sqr_rs + ec0_b2 * rs + ec0_b3 * rs .* sqr_rs + ec0_b4 * rs .* rs);
ec0_q1p = ec0_aa * (ec0_b1 * rsm1_2 + 2.0 * ec0_b2 + 3.0 * ec0_b3 * sqr_rs + 4.0 * ec0_b4 * rs);
ec0_den = 1.0./(ec0_q1 .* ec0_q1 + ec0_q1);
ec0_log = -log(ec0_q1 .* ec0_q1 .* ec0_den);
ecrs0 = ec0_q0 .* ec0_log;
decrs0_drs = -2.0 * ec0_aa * ec0_a1 * ec0_log - ec0_q0 .* ec0_q1p .* ec0_den;

%        Add LSD correlation functional to GGA exchange functional
ec = ecrs0;
v1c = ecrs0 - (rs/3.0) .* decrs0_drs;

%        -----------------------------------------------------------------------------
%        Eventually add the GGA correlation part of the PBE functional
%        Note : the computation of the potential in the spin-unpolarized
%        case could be optimized much further. Other optimizations are left to do.

%        From ec to bb
bb = ecrs0 * gamma_inv;
dbb_drs = decrs0_drs * gamma_inv;

%        From bb to cc
exp_pbe = exp(-bb);
cc = 1.0./(exp_pbe - 1.0);
dcc_dbb = cc .* cc .* exp_pbe;
dcc_drs = dcc_dbb .* dbb_drs;

%        From cc to aa
coeff_aa = beta * gamma_inv;
aa = coeff_aa * cc;
daa_drs = coeff_aa * dcc_drs;

%        Introduce tt : do not assume that the spin-dependent gradients are collinear
dtt_dg = 2.0 * rhotot_inv .* rhotot_inv .* rhom1_3 * coeff_tt;
%        Note that tt is (the t variable of PBE divided by phi) squared
tt = 0.5 * sigma .* dtt_dg;

%        Get xx from aa and tt
xx = aa .* tt;
dxx_drs = daa_drs .* tt;
dxx_dtt = aa;

%        From xx to pade
pade_den = 1.0./(1.0 + xx .* (1.0 + xx));
pade = (1.0 + xx) .* pade_den;
dpade_dxx = -xx .* (2.0 + xx) .* (pade_den.^2);
dpade_drs = dpade_dxx .* dxx_drs;
dpade_dtt = dpade_dxx .* dxx_dtt;

%        From pade to qq
qq = tt .* pade;
dqq_drs = tt .* dpade_drs;
dqq_dtt = pade + tt .* dpade_dtt;

%        From qq to rr
arg_rr = 1.0 + beta * gamma_inv * qq;
div_rr = 1.0./arg_rr;
rr = gamma * log(arg_rr);
drr_dqq = beta * div_rr;
drr_drs = drr_dqq .* dqq_drs;
drr_dtt = drr_dqq .* dqq_dtt;

%        The GGA correlation energy is added
ec = ec + rr;

%        From hh to the derivative of the energy wrt the density
drhohh_drho = rr - third * rs .* drr_drs - (7.0/3.0) * tt .* drr_dtt; %- zeta * dhh_dzeta 
v1c = v1c + drhohh_drho;

%        From hh to the derivative of the energy wrt to the gradient of the
%        density, divided by the gradient of the density
%        (The v3.3 definition includes the division by the norm of the gradient)

v2c = rho .* dtt_dg .* drr_dtt;
end

function [sxsr,v1xsr,v2xsr] = pbexsr(rho,grho,omega)
% pbe short-ranged exchange 

% constants 
us = 0.161620459673995492;
ax = -0.738558766382022406;
% um = 0.2195149727645171;        % mu
% uk = 0.804;                     % kappa
% ul = um/uk;
f1 = -1.10783814957303361;
alpha = 2/3;

fx = zeros(size(rho));
d1x = zeros(size(rho));
d2x = zeros(size(rho));

rs = rho.^(1/3);
vx = (4/3).*f1*alpha*rs;

aa = grho;

rr = 1./(rho.*rs);
ex = ax./rr;
s2 = aa.*rr.*rr*us*us;

s = sqrt(s2);
indx = s > 8.3;
s(indx) = 8.572844 - 18.796223/s2(indx);

for i = 1:size(rho,1)
    for j = 1:size(rho,2)
        [fx(i,j),d1x(i,j),d2x(i,j)] = wpbe_analy_erfc_approx_grad(rho(i,j),s(i,j),omega);
    end
end

sxsr  = ex.*fx;
dsdn  = -4/3*s./rho;
v1xsr = vx.*fx + (dsdn.*d2x+d1x).*ex;
dsdg  = us*rr;
v2xsr = ex./sqrt(aa).*dsdg.*d2x;

end

function [Fx_wpbe,d1rfx,d1sfx] = wpbe_analy_erfc_approx_grad(rho,s,omega)
% Taken from Quantum-Espresso
r36 = 36;
r64 = 64;
r81 = 81;
r256 = 256;
r384 = 384;

r27 = 27;
r128 = 128;
r144 = 144;
r288 = 288;
r324 = 324;
r729  = 729;

r20 = 20;
r32 = 32;
r243 = 243;
r2187 = 2187;
r6561 = 6561;
r40 = 40;

r12 = 12;
r25 = 25;
r30 = 30;
r54 = 54;
r75 = 75;
r105 = 105;
r135 = 135;
r1215 = 1215;
r15309  = 15309;
  
% General constants
f12    = 0.5;
f13    = 1/3;
f14    = 0.25;
f32    = 1.5;
f34    = 0.75;
f94    = 2.25;
f98    = 1.125;
f1516  = 15/16;
pi2    = pi*pi;
srpi   = sqrt(pi);

% Constants from fit
ea1 = -1.128223946706117;
ea2 = 1.452736265762971;
ea3 = -1.243162299390327;
ea4 = 0.971824836115601;
ea5 = -0.568861079687373;
ea6 = 0.246880514820192;
ea7 = -0.065032363850763;
ea8 = 0.008401793031216;
eb1 = 1.455915450052607;

% Constants for PBE hole
A      =  1.0161144;
B      = -3.7170836e-1;
C      = -7.7215461e-2;
D      =  5.7786348e-1;
E      = -5.1955731e-2;

% Constants for fit of H(s) (PBE)
Ha1    = 9.79681e-3;
Ha2    = 4.10834e-2;
Ha3    = 1.87440e-1;
Ha4    = 1.20824e-3;
Ha5    = 3.47188e-2;

% Constants for F(H) (PBE)
Fc1    = 6.4753871;
Fc2    = 4.7965830e-1;

% Constants for polynomial expansion for EG for small s
EGa1   = -2.628417880e-2;
EGa2   = -7.117647788e-2;
EGa3   =  8.534541323e-2;

% Constants for large x expansion of exp(x)*ei(-x)
expei1 = 4.03640;
expei2 = 1.15198;
expei3 = 5.03627;
expei4 = 4.19160;

% Cutoff criterion below which to use polynomial expansion
EGscut     = 8.0e-2;
wcutoff    = 1.4e1;
expfcutoff = 7.0e2;

xkf    = (3*pi2*rho) ^ f13;

A2 = A*A;
A3 = A2*A;
A12 = sqrt(A);
A32 = A12*A;
A52 = A32*A;

w     = omega / xkf;
w2    = w * w;
w3    = w2 * w;
w4    = w2 * w2;
w5    = w3 * w2;
w6    = w5 * w;
w7    = w6 * w;
w8    = w7 * w;
d1rw  = -(1/(3*rho))*w;
X      = - 8/9;

s2     = s*s;
s3     = s2*s;
s4     = s2*s2;
s5     = s4*s;
s6     = s5*s;

% Calculate wPBE enhancement factor;

Hnum    = Ha1*s2 + Ha2*s4;
Hden    = 1 + Ha3*s4 + Ha4*s5 + Ha5*s6;

H       = Hnum/Hden;

d1sHnum = 2*Ha1*s + 4*Ha2*s3;
d1sHden = 4*Ha3*s3 + 5*Ha4*s4 + 6*Ha5*s5;

d1sH    = (Hden*d1sHnum - Hnum*d1sHden) / (Hden*Hden);

F      = Fc1*H + Fc2;
d1sF   = Fc1*d1sH;

% Change exp1nt of Gaussian if we're using the simple approx.;

if w > wcutoff 
    eb1 = 2.0d0;
end

% Calculate helper variables (should be moved later on...);

Hsbw = s2*H + eb1*w2;
Hsbw2 = Hsbw*Hsbw;
Hsbw3 = Hsbw2*Hsbw;
Hsbw4 = Hsbw3*Hsbw;
Hsbw12 = sqrt(Hsbw);
Hsbw32 = Hsbw12*Hsbw;
Hsbw52 = Hsbw32*Hsbw;
Hsbw72 = Hsbw52*Hsbw;

d1sHsbw  = d1sH*s2 + 2*s*H;
d1rHsbw  = 2*eb1*d1rw*w;

DHsbw = D + s2*H + eb1*w2;
DHsbw2 = DHsbw*DHsbw;
DHsbw3 = DHsbw2*DHsbw;
DHsbw4 = DHsbw3*DHsbw;
DHsbw5 = DHsbw4*DHsbw;
DHsbw12 = sqrt(DHsbw);
DHsbw32 = DHsbw12*DHsbw;
DHsbw52 = DHsbw32*DHsbw;
DHsbw72 = DHsbw52*DHsbw;
DHsbw92 = DHsbw72*DHsbw;

HsbwA94   = f94 * Hsbw / A;
HsbwA942  = HsbwA94*HsbwA94;
HsbwA943  = HsbwA942*HsbwA94;
HsbwA945  = HsbwA943*HsbwA942;
HsbwA9412 = sqrt(HsbwA94);

DHs    = D + s2*H;
DHs2   = DHs*DHs;
DHs3   = DHs2*DHs;
DHs4   = DHs3*DHs;
DHs72  = DHs3*sqrt(DHs);
DHs92  = DHs72*DHs;

d1sDHs = 2*s*H + s2*d1sH;

DHsw   = DHs + w2;
DHsw2  = DHsw*DHsw;
DHsw52 = sqrt(DHsw)*DHsw2;
DHsw72 = DHsw52*DHsw;

d1rDHsw = 2*d1rw*w;

if s > EGscut
    G_a    = srpi*(15*E+6*C*(1+F*s2)*DHs+4*B*(DHs2)+8*A*(DHs3))*(1/(16*DHs72))...
             - f34*pi*sqrt(A) * exp(f94*H*s2/A)*(1 - erf(f32*s*sqrt(H/A)));
    d1sG_a = (1/r32)*srpi * ((r36*(2*H + d1sH*s) / (A12*sqrt(H/A)))+ (1/DHs92) ...
              * (-8*A*d1sDHs*DHs3 - r105*d1sDHs*E-r30*C*d1sDHs*DHs*(1+s2*F)...
              +r12*DHs2*(-B*d1sDHs + C*s*(d1sF*s + 2*F)))-((r54*exp(f94*H*s2/A)...
              *srpi*s*(2*H+d1sH*s)*erfc(f32*sqrt(H/A)*s))/ A12));
    G_b    = (f1516 * srpi * s2) / DHs72;
    d1sG_b = (15*srpi*s*(4*DHs - 7*d1sDHs*s))/(r32*DHs92);
    EG     = - (f34*pi + G_a) / G_b;
    d1sEG  = (-4*d1sG_a*G_b + d1sG_b*(4*G_a + 3*pi))/(4*G_b*G_b);

else
    EG    = EGa1 + EGa2*s2 + EGa3*s4;
    d1sEG = 2*EGa2*s + 4*EGa3*s3;
    
end

% Calculate the terms needed in any case
      
term2 = (DHs2*B + DHs*C + 2*E + DHs*s2*C*F + 2*s2*EG)/(2*DHs3);

d1sterm2 = (-6*d1sDHs*(EG*s2 + E) + DHs2 * (-d1sDHs*B + s*C*(d1sF*s + 2*F)) ...
           + 2*DHs * (2*EG*s - d1sDHs*C + s2 * (d1sEG - d1sDHs*C*F))) / (2*DHs4);

term3 = - w  * (4*DHsw2*B + 6*DHsw*C + 15*E + 6*DHsw*s2*C*F + 15*s2*EG) / (8*DHs*DHsw52);

d1sterm3 = w * (2*d1sDHs*DHsw * (4*DHsw2*B + 6*DHsw*C + 15*E + 3*s2*(5*EG + 2*DHsw*C*F)) ...
           + DHs * (r75*d1sDHs*(EG*s2 + E) + 4*DHsw2*(d1sDHs*B - 3*s*C*(d1sF*s + 2*F))   ...
           - 6*DHsw*(-3*d1sDHs*C + s*(10*EG + 5*d1sEG*s - 3*d1sDHs*s*C*F))))    ...
           / (16*DHs2*DHsw72);

d1rterm3 = (-2*d1rw*DHsw * (4*DHsw2*B + 6*DHsw*C + 15*E + 3*s2*(5*EG + 2*DHsw*C*F)) ...
           + w * d1rDHsw * (r75*(EG*s2 + E) + 2*DHsw*(2*DHsw*B + 9*C + 9*s2*C*F))) ...
           / (16*DHs*DHsw72);

term4 = - w3 * (DHsw*C + 5*E + DHsw*s2*C*F + 5*s2*EG) / (2*DHs2*DHsw52);

d1sterm4 = (w3 * (4*d1sDHs*DHsw * (DHsw*C + 5*E + s2 * (5*EG + DHsw*C*F))...
           + DHs * (r25*d1sDHs*(EG*s2 + E) - 2*DHsw2*s*C*(d1sF*s + 2*F) ...
           + DHsw * (3*d1sDHs*C + s*(-r20*EG - 10*d1sEG*s + 3*d1sDHs*s*C*F)))))...
           / (4*DHs3*DHsw72);

d1rterm4 = (w2 * (-6*d1rw*DHsw * (DHsw*C + 5*E + s2 * (5*EG + DHsw*C*F))...
          + w * d1rDHsw * (r25*(EG*s2 + E) + 3*DHsw*C*(1 + s2*F))))  ...
          / (4*DHs2*DHsw72);

term5 = - w5 * (E + s2*EG) / (DHs3*DHsw52);

d1sterm5 = (w5 * (6*d1sDHs*DHsw*(EG*s2 + E) + DHs * (-2*DHsw*s * ...
           (2*EG + d1sEG*s) + 5*d1sDHs * (EG*s2 + E)))) / (2*DHs4*DHsw72);

d1rterm5 = (w4 * 5*(EG*s2 + E) * (-2*d1rw*DHsw + d1rDHsw * w)) / (2*DHs3*DHsw72);



if s > 0.0 || w > 0.0
    t10    = (f12)*A*log(Hsbw / DHsbw);
    t10d1  = f12*A*(1/Hsbw - 1/DHsbw);
    d1st10 = d1sHsbw*t10d1;
    d1rt10 = d1rHsbw*t10d1;
end
  
% Calculate exp(x)*f(x) depending on size of x
if HsbwA94 < expfcutoff 
    piexperf = pi*exp(HsbwA94)*erfc(HsbwA9412);
    expei    = exp(HsbwA94)*(-expint(HsbwA94));
else
    piexperf = pi*(1/(srpi*HsbwA9412) - 1/(2*sqrt(pi*HsbwA943)) + 3/(4*sqrt(pi*HsbwA945)));
    
    expei  = - (1/HsbwA94) * (HsbwA942 + expei1*HsbwA94 + expei2)/(HsbwA942 + expei3*HsbwA94 + expei4);
end

% Calculate the derivatives (based on the orig. expression)
piexperfd1  = - (3*srpi*sqrt(Hsbw/A))/(2*Hsbw) + (9*piexperf)/(4*A);
d1spiexperf = d1sHsbw*piexperfd1;
d1rpiexperf = d1rHsbw*piexperfd1;

expeid1  = f14*(4/Hsbw + (9*expei)/A);
d1sexpei = d1sHsbw*expeid1;
d1rexpei = d1rHsbw*expeid1;



if w == 0
    % Fall back to original expression for the PBE hole
    t1 = -f12*A*expei;
    d1st1 = -f12*A*d1sexpei;
    d1rt1 = -f12*A*d1rexpei;
    
    if s > 0.0
        term1    = t1 + t10;
        d1sterm1 = d1st1 + d1st10;
        d1rterm1 = d1rt1 + d1rt10;
        Fx_wpbe = X * (term1 + term2);
        d1sfx = X * (d1sterm1 + d1sterm2);
        d1rfx = X * d1rterm1;
    else
        Fx_wpbe = 1.0d0;
        % TODO    This is checked to be true for term1
        %         How about the other terms???
        d1sfx   = 0.0d0;
        d1rfx   = 0.0d0;
    end
    
elseif w > wcutoff
    % Use simple Gaussian approximation for large w
    term1 = -f12*A*(expei+log(DHsbw)-log(Hsbw));
    term1d1  = - A/(2*DHsbw) - f98*expei;
    d1sterm1 = d1sHsbw*term1d1;
    d1rterm1 = d1rHsbw*term1d1;

    Fx_wpbe = X * (term1 + term2 + term3 + term4 + term5);
    d1sfx = X * (d1sterm1 + d1sterm2 + d1sterm3 + d1sterm4 + d1sterm5);
    d1rfx = X * (d1rterm1 + d1rterm3 + d1rterm4 + d1rterm5);

else 
    % For everything else, use the full blown expression
    % First, we calculate the polynomials for the first term
    np1    = -f32*ea1*A12*w + r27*ea3*w3/(8*A12) - r243*ea5*w5/(r32*A32) + r2187*ea7*w7/(r128*A52);
    d1rnp1 = - f32*ea1*d1rw*A12 + (r81*ea3*d1rw*w2)/(8*A12) - ...
             (r1215*ea5*d1rw*w4)/(r32*A32) + (r15309*ea7*d1rw*w6)/(r128*A52);
    np2 = -A + f94*ea2*w2 - r81*ea4*w4/(16*A) + r729*ea6*w6/(r64*A2) - r6561*ea8*w8/(r256*A3);
    d1rnp2 =   f12*(9*ea2*d1rw*w) - (r81*ea4*d1rw*w3)/(4*A)    ...
             + (r2187*ea6*d1rw*w5)/(r32*A2) - (r6561*ea8*d1rw*w7)/(r32*A3);

    % The first term is
    t1    = f12*(np1*piexperf + np2*expei);
    d1st1 = f12*(d1spiexperf*np1 + d1sexpei*np2);
    d1rt1 = f12*(d1rnp2*expei + d1rpiexperf*np1 + d1rexpei*np2 + d1rnp1*piexperf);

    % The factors for the main polynomoal in w and their derivatives
    f2    = (f12)*ea1*srpi*A / DHsbw12;
    f2d1  = - ea1*srpi*A / (4*DHsbw32);
    d1sf2 = d1sHsbw*f2d1;
    d1rf2 = d1rHsbw*f2d1;

    f3    = (f12)*ea2*A / DHsbw;
    f3d1  = - ea2*A / (2*DHsbw2);
    d1sf3 = d1sHsbw*f3d1;
    d1rf3 = d1rHsbw*f3d1;

    f4    =  ea3*srpi*(-f98 / Hsbw12 + f14*A / DHsbw32);
    f4d1  = ea3*srpi*((9/(16*Hsbw32))- (3*A/(8*DHsbw52)));
    d1sf4 = d1sHsbw*f4d1;
    d1rf4 = d1rHsbw*f4d1;

    f5    = ea4*(1/r128) * (-r144*(1/Hsbw) + r64*(1/DHsbw2)*A);
    f5d1  = ea4*((f98/Hsbw2)-(A/DHsbw3));
    d1sf5 = d1sHsbw*f5d1;
    d1rf5 = d1rHsbw*f5d1;

    f6    = ea5*(3*srpi*(3*DHsbw52*(9*Hsbw-2*A) + 4*Hsbw32*A2)) / (r32*DHsbw52*Hsbw32*A);
    f6d1  = ea5*srpi*((r27/(r32*Hsbw52))-(r81/(r64*Hsbw32*A))-((15*A)/(16*DHsbw72)));
    d1sf6 = d1sHsbw*f6d1;
    d1rf6 = d1rHsbw*f6d1;

    f7    = ea6*(((r32*A)/DHsbw3 + (-r36 + (r81*s2*H)/A)/Hsbw2)) / r32;
    d1sf7 = ea6*(3*(r27*d1sH*DHsbw4*Hsbw*s2 + 8*d1sHsbw*A*(3*DHsbw4 - 4*Hsbw3*A) + ...
            r54*DHsbw4*s*(Hsbw - d1sHsbw*s)*H))/(r32*DHsbw4*Hsbw3*A);
    d1rf7 = ea6*d1rHsbw*((f94/Hsbw3)-((3*A)/DHsbw4)-((r81*s2*H)/(16*Hsbw3*A)));

    f8    = ea7*(-3*srpi*(-r40*Hsbw52*A3+9*DHsbw72*(r27*Hsbw2-6*Hsbw*A+4*A2))) / (r128 * DHsbw72*Hsbw52*A2);
    f8d1  = ea7*srpi*((r135/(r64*Hsbw72)) + (r729/(r256*Hsbw32*A2))  ...
            -(r243/(r128*Hsbw52*A))-((r105*A)/(r32*DHsbw92)));
    d1sf8 = d1sHsbw*f8d1;
    d1rf8 = d1rHsbw*f8d1;

    f9    = (r324*ea6*eb1*DHsbw4*Hsbw*A + ea8*(r384*Hsbw3*A3 + DHsbw4*(-r729*Hsbw2 ...
            + r324*Hsbw*A - r288*A2))) / (r128*DHsbw4*Hsbw3*A2);
    f9d1  = -((r81*ea6*eb1)/(16*Hsbw3*A)) + ea8*((r27/(4*Hsbw4))+(r729/(r128*Hsbw2*A2)) ...
            -(r81/(16*Hsbw3*A))-((r12*A/DHsbw5)));
    d1sf9 = d1sHsbw*f9d1;
    d1rf9 = d1rHsbw*f9d1;

    t2t9    = f2*w  + f3*w2 + f4*w3 + f5*w4 + f6*w5 + f7*w6 + f8*w7 + f9*w8;
    d1st2t9 = d1sf2*w + d1sf3*w2 + d1sf4*w3 + d1sf5*w4 ...
              + d1sf6*w5 + d1sf7*w6 + d1sf8*w7 + d1sf9*w8;
    d1rt2t9 = d1rw*f2 + d1rf2*w + 2*d1rw*f3*w + d1rf3*w2 + 3*d1rw*f4*w2 ...
            + d1rf4*w3 + 4*d1rw*f5*w3 + d1rf5*w4 + 5*d1rw*f6*w4 ...
            + d1rf6*w5 + 6*d1rw*f7*w5 + d1rf7*w6 + 7*d1rw*f8*w6 ...
            + d1rf8*w7 + 8*d1rw*f9*w7 + d1rf9*w8;

    % The final value of term1 for 0 < omega < wcutoff is:
    term1 = t1 + t2t9 + t10;
    d1sterm1 = d1st1 + d1st2t9 + d1st10;
    d1rterm1 = d1rt1 + d1rt2t9 + d1rt10;
    
    % The final value for the enhancement factor and its derivatives is:
    Fx_wpbe = X * (term1 + term2 + term3 + term4 + term5);
    d1sfx = X * (d1sterm1 + d1sterm2 + d1sterm3 + d1sterm4 + d1sterm5);
    d1rfx = X * (d1rterm1 + d1rterm3 + d1rterm4 + d1rterm5);
end

end
