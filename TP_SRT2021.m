%% Initialisation
clear all;
close all;
clc

k = 8*188; % nombre de bits transmis
nb = 2; % nombre de bits par symbole
Ds = 250e3; % debit symbole
Fse = 4; % différence entre le débit symbole et la fréquence d'échantillonage
Fe = Fse * Ds; % fréquence d'échantillonage
span = 8; % largeur du filtre

g = comm.RaisedCosineTransmitFilter('RolloffFactor',0.35,'FilterSpanInSymbols',span,'OutputSamplesPerSymbol',Fse); % filtre de mise en forme
ga = comm.RaisedCosineReceiveFilter('RolloffFactor',0.35,'FilterSpanInSymbols',span,'InputSamplesPerSymbol',Fse,'DecimationFactor',Fse); % filtre adapte

calculError = comm.ErrorRate('ReceiveDelay', 2*span, 'ComputationDelay', 2);



%% Emetteur

bk = double ((randn(1,k)>0).'); % generation des bits
qpsk_mod = comm.QPSKModulator('PhaseOffset',0,'BitInput',1);
an = qpsk_mod.step(bk);

dn = zeros(size(an));
dn(1) = an(1);
for n = 2:size(an)
   dn(n) = an(n)*dn(n-1);
end

dn_se = dn;

sl = step(g,dn_se);

%% Canal

sl_in = sl;

%% Recepteur

rn_se = step(ga,sl_in);
rn = rn_se; %downsample(rn_se,Fse,4); % pas sur

vn = zeros(size(rn_se));
vn(1) = rn_se(1);
for n = 2:size(an)
   vn(n) = rn_se(n) * conj(rn_se(n-1));
end

qpsk_demod = comm.QPSKDemodulator('PhaseOffset',0,'BitOutput',1);
anr = qpsk_demod.step(vn);
bkr = anr;

error_cnt = calculError.step( bk, bkr );
TEB = error_cnt(1)

%% Figures

%constellation(qpsk_mod);
%constellation(qpsk_demod);
%figure
%plot(dn, '*');
figure;
plot(vn, '*');
























