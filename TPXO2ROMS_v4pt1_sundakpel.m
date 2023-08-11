%Invoke function TPXO2ROMS_v4pt1.

%lengthSim is the approximate anticipated length of model run (days) (nodal corrections are computed at mid-point).
%fnGrid is your ROMS grid file.
%TPXO2ROMS will write to fnOut.
%The available choices for ROMSnames are:  'MM' 'MF' 'Q1' 'O1' 'P1' 'K1' 'N2' 'M2' 'S2' 'K2' 'MN4' 'M4' 'MS4'
%For example, use ROMSnames={'M2' 'S2' 'MM'}; to use M2, S2, and MM only
%MM, MF, MN4, and MS4 (the 1/6 degree resolution harmonics) are not well-tested (I don't trust them).

clear
clc
clear all

%****************Modify to suit***************

%% cara 1 : 
% tide start diisi dengan jumlah hari selisis tanggal awal simulasi
% dengan tanggal awal tahun simulasi; lengthsim diisidengan jumlah hari
% simulasil; t0 diisi dengan tanggal awal tahun simulasi

% TIDE_START=273;  %ROMS daynumber in ROMS *.in file (TIME_REF=-2) selisis tanggal awal simulasi dengan t0
% lengthSim=7;  %estimated length of model run in days (not exact)
% t0=datenum(2022, 10, 1, 0, 0, 0); % tanggal awal tahun waktu simulasi

%% cara 2 :
% tidestart diisi dengan jumlah hari selisis tanggal awal simulasi dengan
% tanggal awal julian day (1968-05-23); lengthsim diisidengan jumlah hari
% simulasi; t0 diisi dengan tanggal awal tahun simulasi (julian day
% (1968-05-23))

%TIDE_START=datenum(2022,10,1)-datenum(1968,05,23);  % datenum(2022,10,1)-datenum(1968,05,23)
%lengthSim=7;  %estimated length of model run in days (not exact)
%t0=datenum(1968, 05, 23, 0, 0, 0); % tanggal awal tahun waktu simulasi

%% cara 3 :
% tidestart diisi dengan jumlah hari selisis tanggal awal simulasi dengan
% tanggal ref di script roms_master_climatology_coawst_mw (1858,11,17); 
% lengthsim diisidengan jumlah hari simulasi; 
% t0 diisi dengan tanggal awal tahun simulasi (1858,11,17)

TIDE_START=datenum(2022,12,29)-datenum(1858,11,17);  % datenum(2022,10,1)-datenum(1968,05,23)
lengthSim=3;  %estimated length of model run in days (not exact)
t0=datenum(1858, 11, 17, 0, 0, 0); % tanggal awal tahun waktu simulasi

%% set input dan output file
fnGrid='C:\Arief\ROMS\03_GridBuilder\01_Selat_Sunda\Kopel_Sunda3\grid_sundaD_9kmv3_kopel_batnas.nc'; 
%fnOut=['C:\Arief\ROMS\08_Tidal\TPXO8\02_Sunda\tides_sundaD.nc'];
fnOut=['C:\Arief\ROMS\08_Tidal\TPXO8\03_SundaKopel\tides_sundakopelv3.nc'];
ROMSnames={'Q1' 'O1' 'P1' 'K1' 'N2' 'M2' 'S2' 'K2'};


%% ********************************************

%t0=TIDE_START+datenum(1968,5,23);  %Convert TIDE_START to Matlab daynumber

disp(['Output filename: ' fnOut])

%Check file
if exist(fnOut,'file');
   error('File fnOut exists - must delete it or modify fnOut filename.')
end
%Invoke
TPXO2ROMS_v4pt1(t0,ROMSnames,fnGrid,fnOut,lengthSim)
