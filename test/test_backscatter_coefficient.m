% Sat 21 Jul 15:13:47 CEST 2018

f   = 1.2e6;
d_mm= 0.1;

% ssc=2.650*100*1e-6;
% 100mg/l = 100*1e-3kg/m^3
% ssc  = 100*1e-3
ssc    = 100*1e-3

bs   = ssc2backscatter(ssc,d_mm,f)
ssc_ = bs/backscatter_coefficient_2(d_mm,f)

'test'
ssc/ssc_

Sv = 10*log10(bs)

