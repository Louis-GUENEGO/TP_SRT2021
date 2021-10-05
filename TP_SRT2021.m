%% Initialisation
clear all;
close all;
clc

k = 8*188; % nombre de bits transmis
nb = 2; % nombre de bits par symbole
Ds = 250e3; % debit symbole
Fse = 4;
Fe = Fse * Ds;

g = comm.RaisedCosineTransmitFilter('RolloffFactor',0.35,'FilterSpanInSymbols',8,'OutputSamplesPerSymbol',Fse); % filtre de mise en forme
ga = comm.RaisedCosineReceiveFilter('RolloffFactor',0.35,'FilterSpanInSymbols',8,'InputSamplesPerSymbol',Fse,'DecimationFactor',Fse); % filtre adapte

%% Emetteur

bk = (randn(1,k)>0).'; % generation des bits
qpsk_mod = comm.QPSKModulator('PhaseOffset',0,'BitInput',1);
an = qpsk_mod.step(bk);

dn = an.*dn(n-1); %boucle for ou  cumprod?

dn_se = dn;

sl = g(dn_se);

%% Canal

sl_in = sl;

%% Recepteur

rn_se = ga(sl_in);
rn = rn_se;%downsample(rn_se,Fse,4); % pas sur
vn = rn;

qpsk_demod = comm.QPSKDemodulator('PhaseOffset',0,'BitOutput',1);
anr = qpsk_demod.step(vn);
bkr = anr;

offset = 16;
mat1 = bkr(offset+1:k,1);
mat2 = double(bk(1:k-offset,1));

error_cnt = biterr(mat1, mat2);
TEB = error_cnt/(k-offset);
%comm.ErrorRate();


%% Figures

%constellation(qpsk_mod);
%constellation(qpsk_demod);
%plot(an, '*');
%figure;
%plot(dn, '*');
























