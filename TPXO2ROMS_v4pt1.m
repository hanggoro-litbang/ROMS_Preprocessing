function TPXO2ROMS_v4pt1(t0,ROMSnames,fnGrid,fnOut,ndays)

% function TPXO2ROMS_v4pt1(t0,ROMSnames,fnGrid,fnOut,ndays)
% t0: start time
% ROMSnames: ROMS components to calculate
% typical components: MM  MF  Q1  O1  P1  K1  N2  M2  S2  K2  MN4 M4  MS4
% fnGrid: ROMS grid to calculate tidal components on
% fnOut: name of ROMS tide file
% ndays: Approximate length of run in days (default 365)

%Prepares a tidal forcing file for ROMS from TPXO tidal model (OSU) version 8.
%This version requires Matlab R2012a or later because it takes advantage of
%TriScatteredInterp being able to use complex values from R2012a on.
%If run on R2013 on it will use scatteredInterpolant/griddedInterpolant
%instead of TriScatteredInterp which will be discontinued at some point.

%Uses data from: http://volkov.oce.orst.edu/tides/tpxo8_atlas.html
%Download the netcdf files and put them in a folder on your Matlab path.

%Notes for version 4.1
%20 May 2015
%Minor improvements to code, e.g. interpolation choice.

%Notes for version 4.0
%MM, MF, MS4, and MN4 (the 1/6th degree harmonics) don't look quite right. I
%would avoid them until resolved.

%Notes for version 3.1
%20 October 2013: bug fixes (saves velocity, not
%transport; added mod(lon,360) for cases where lon ran from
%-180 to 180 instead of from 0 to 360, but note that the interpolation may
%create a small distortion close to the Greenwich meridian (between 359-360
%degrees east).

%Notes for version 3.0
%4 October 2013: has Charles James' improvements (thanks Charles!).
%Uses Matlab native calls instead of snctools - not sure I consider
%this an improvement, but it means users don't have to install snctools;
%it also means that arrays are now ordered as column, row, Nth harmonic);
%improved interpolation; faster; cleaner.

%Notes for version 2.0
%14 August 2012: fixed for case of
%extracting a single component; got rid of PRESERVE_FVD and related
%comments; modified interpTPXO.m.

%Created by J. Luick
%www.austides.com 23October2013
%**************************************************************************

varcheck('ndays',365);    %Sets default ndays=365 if none specified.
varcheck('ROMSnames',{'MM','MF','Q1','O1','P1','K1','N2','M2','S2','K2','MN4','M4','MS4'});
varcheck('fnOut','Ocean_tide.nc');

lengthSim=ndays;     %Approximate anticipated length of model run (days) (for f & u)
ROMStitle=['ROMS TPXO data for ' datestr(t0,0)];

%ROMS grid info
lonR=ncread(fnGrid,'lon_rho');
lonR=mod(lonR,360); %This doesn't solve the problem of grids that span longitude=0;
latR=ncread(fnGrid,'lat_rho');
maskR=ncread(fnGrid,'mask_rho');
maskP=ncread(fnGrid,'mask_psi');
[L,M]=size(maskP);

%TPXO filenames
TPXOfile.grid.atlas6='grid_tpxo8_atlas6.nc';
TPXOfile.grid.atlas30='grid_tpxo8atlas_30.nc';

TPXOfile.elev.k1='hf.k1_tpxo8_atlas_30c.nc';
TPXOfile.elev.k2='hf.k2_tpxo8_atlas_30c.nc';
TPXOfile.elev.m2='hf.m2_tpxo8_atlas_30c.nc';
TPXOfile.elev.m4='hf.m4_tpxo8_atlas_30c.nc';
TPXOfile.elev.mf='hf.mf_tpxo8_atlas_6.nc';
TPXOfile.elev.mm='hf.mm_tpxo8_atlas_6.nc';
TPXOfile.elev.mn4='hf.mn4_tpxo8_atlas_6.nc';
TPXOfile.elev.ms4='hf.ms4_tpxo8_atlas_6.nc';
TPXOfile.elev.n2='hf.n2_tpxo8_atlas_30c.nc';
TPXOfile.elev.o1='hf.o1_tpxo8_atlas_30c.nc';
TPXOfile.elev.p1='hf.p1_tpxo8_atlas_30c.nc';
TPXOfile.elev.q1='hf.q1_tpxo8_atlas_30c.nc';
TPXOfile.elev.s2='hf.s2_tpxo8_atlas_30c.nc';

TPXOfile.vel.k1='uv.k1_tpxo8_atlas_30c.nc';
TPXOfile.vel.k2='uv.k2_tpxo8_atlas_30c.nc';
TPXOfile.vel.m2='uv.m2_tpxo8_atlas_30c.nc';
TPXOfile.vel.m4='uv.m4_tpxo8_atlas_30c.nc';
TPXOfile.vel.mf='uv.mf_tpxo8_atlas_6.nc';
TPXOfile.vel.mm='uv.mm_tpxo8_atlas_6.nc';
TPXOfile.vel.mn4='uv.mn4_tpxo8_atlas_6.nc';
TPXOfile.vel.ms4='uv.ms4_tpxo8_atlas_6.nc';
TPXOfile.vel.n2='uv.n2_tpxo8_atlas_30c.nc';
TPXOfile.vel.o1='uv.o1_tpxo8_atlas_30c.nc';
TPXOfile.vel.p1='uv.p1_tpxo8_atlas_30c.nc';
TPXOfile.vel.q1='uv.q1_tpxo8_atlas_30c.nc';
TPXOfile.vel.s2='uv.s2_tpxo8_atlas_30c.nc';

%TPXO data
disp('Reading TPXO')
TPXO=readTPXOdata(TPXOfile,ROMSnames,lonR,latR);

%Magic numbers for TPXO (Solar Doodson Numbers, Ref. Phase, Speed degrees/hour)
TPXOharmonics.MM= [0   1   0  -1   0   0    0    0.5443747];
TPXOharmonics.MF= [0   2   0   0   0   0    0    1.0980331];
TPXOharmonics.Q1= [1  -3   1   1   0   0  270   13.3986607];
TPXOharmonics.O1= [1  -2   1   0   0   0  270   13.9430351];
TPXOharmonics.P1= [1   0  -1   0   0   0  270   14.9589310];
TPXOharmonics.K1= [1   0   1   0   0   0   90   15.0410690];
TPXOharmonics.N2= [2  -3   2   1   0   0    0   28.4397297];
TPXOharmonics.M2= [2  -2   2   0   0   0    0   28.9841042];
TPXOharmonics.S2= [2   0   0   0   0   0    0   30.0000000];
TPXOharmonics.K2= [2   0   2   0   0   0    0   30.0821381];
TPXOharmonics.MN4=[4  -5   4   1   0   0    0   57.4238319];
TPXOharmonics.M4= [4  -4   4   0   0   0    0   57.9682083];
TPXOharmonics.MS4=[4  -2   2   0   0   0    0   58.9841042];

%Magic numbers the way ROMS sees them
Nharmonics=length(ROMSnames);
ROMSperiods=zeros(1,Nharmonics);
for n=1:Nharmonics
    ROMSharmonics.(ROMSnames{n})=TPXOharmonics.(ROMSnames{n});
    ROMSperiods(n)=360./TPXOharmonics.(ROMSnames{n})(end);
    Vdeg.(ROMSnames{n})=Vphase(t0,ROMSharmonics.(ROMSnames{n}));
end

% V,u,f (reference phase and nodal corrections)
[fFac,uFac]=TPXOnodalfactors(t0+lengthSim/2,ROMSnames);

% Extract tide info from TPXO and put on rho grid
zamp=zeros(L+1,M+1,Nharmonics);
zpha=zeros(L+1,M+1,Nharmonics);
uamp=zeros(L+1,M+1,Nharmonics);
upha=zeros(L+1,M+1,Nharmonics);
vamp=zeros(L+1,M+1,Nharmonics);
vpha=zeros(L+1,M+1,Nharmonics);
major=zeros(L+1,M+1,Nharmonics);
eccentricity=zeros(L+1,M+1,Nharmonics);
inclination=zeros(L+1,M+1,Nharmonics);
phase=zeros(L+1,M+1,Nharmonics);

count_a6=0;
count_a30=0;
for k=1:Nharmonics
   if strcmp(ROMSnames{k},'MS4') || strcmp(ROMSnames{k},'MN4') || ...
      strcmp(ROMSnames{k},'MF') || strcmp(ROMSnames{k},'MM')
      count_a6=count_a6+1;
      count=count_a6;
      harmonic=TPXO.harmonic_a6.harmonic(count_a6);
   else
      count_a30=count_a30+1;
      count=count_a30;
      harmonic=TPXO.harmonic_a30.harmonic(count_a30);
   end
   disp(['Interpolating ',harmonic,' amplitudes'])
   ei=interpTPXO(TPXO.h,harmonic,count,lonR,latR,maskR);
   zamp(:,:,k)=abs(ei).*fFac.(ROMSnames{k});
   zpha(:,:,k)=mod(-angle(ei)*180/pi-uFac.(ROMSnames{k})-Vdeg.(ROMSnames{k}),360);

   disp(['Interpolating ',harmonic,' u components'])
   ei=interpTPXO(TPXO.U,harmonic,count,lonR,latR,maskR);
   uamp(:,:,k)=abs(ei).*fFac.(ROMSnames{k});
   upha(:,:,k) =mod(-angle(ei)*180/pi-uFac.(ROMSnames{k})-Vdeg.(ROMSnames{k}),360);

   disp(['Interpolating ',harmonic,' v components'])
   ei=interpTPXO(TPXO.V,harmonic,count,lonR,latR,maskR);
   vamp(:,:,k)=abs(ei).*fFac.(ROMSnames{k});
   vpha(:,:,k)=mod(-angle(ei)*180/pi-uFac.(ROMSnames{k})-Vdeg.(ROMSnames{k}),360);

   [maj,ecc,inc,pha]=ap2ep(squeeze(uamp(:,:,k)),squeeze(upha(:,:,k)),...
      squeeze(vamp(:,:,k)),squeeze(vpha(:,:,k)));
   ecc(isnan(ecc))=0; %zero current results in e=NaN. ROMS crashes. (ps circle has e=0)
   major(:,:,k)=maj;
   eccentricity(:,:,k)=ecc;
   inclination(:,:,k)=inc;
   phase(:,:,k)=pha;
end

minor=major.*eccentricity;

%Set up output netcdf forcing file

    %global atts
    %String version of ROMSnames (for global attribute only)
    varname={'tide_period'; 'tide_Ephase'; 'tide_Eamp';...
        'tide_Cmin'; 'tide_Cmax'; 'tide_Cangle'; 'tide_Cphase';...
        'tide_Uamp'; 'tide_Uphase'; 'tide_Vamp'; 'tide_Vphase'};

    nVar=length(varname);

    if ~exist(fnOut,'file');
        ROMSnames_att=cellfun(@(x) [x, ' '],ROMSnames,'uniformoutput',false);
        ROMSnames_att=[ROMSnames_att{:}];

        attvals={{'Tide angular period','hours'};           %tide_period
            {'Tide elevation phase angle','degrees'}; %tide_Ephase
            {'Tide elevation amplitude','meters'};    %tide_Eamp
            {'Tidal current ellipse semi-minor axis','meter second-1'}; %tide_Cmin
            {'Tidal current ellipse semi-major axis','meter second-1'}; %tide_Cmax
            {'Tidal current ellipse inclination angle','degrees between semi-major axis and east'}; %tide_Cangle
            {'Tidal current phase angle','degrees'};                  %tide_Cphase
            {'Tidal current U-component amplitude','meters'}  %tide_Uamp
            {'Tidal current U-component phase','degrees'}        %tide_Uphase
            {'Tidal current V-component amplitude','meters'}   %tide_Vamp
            {'Tidal current V-component phase','degrees'}};      %tide_Vphase

        fileschema.Name='/';
        fileschema.Dimensions(1).Name='xi_rho';
        fileschema.Dimensions(1).Length=L+1;
        fileschema.Dimensions(2).Name='eta_rho';
        fileschema.Dimensions(2).Length=M+1;
        fileschema.Dimensions(3).Name='tide_period';
        fileschema.Dimensions(3).Length=Nharmonics;
        fileschema.Attributes(1).Name='title';
        fileschema.Attributes(1).Value=ROMStitle;
        fileschema.Attributes(2).Name='Creation_date';
        fileschema.Attributes(2).Value=datestr(date,'yyyymmdd');
        fileschema.Attributes(3).Name='grd_file';
        fileschema.Attributes(3).Value=fnGrid;
        fileschema.Attributes(4).Name='type';
        fileschema.Attributes(4).Value='ROMS forcing file from TPXO';
        fileschema.Attributes(5).Name='ini_date_datenumber';
        fileschema.Attributes(5).Value=t0;
        fileschema.Attributes(6).Name='ini_date_mjd';
        fileschema.Attributes(6).Value=t0-datenum(1968,5,23);
        fileschema.Attributes(7).Name='components';
        fileschema.Attributes(7).Value=ROMSnames_att;
        fileschema.Format='64bit';
        for i=1:nVar
            fileschema.Variables(i).Name=varname{i};
            if strcmp(varname{i},'tide_period');
                fileschema.Variables(i).Dimensions(1).Name='tide_period';
                fileschema.Variables(i).Dimensions(1).Length=Nharmonics;
            else
                fileschema.Variables(i).Dimensions(1).Name='xi_rho';
                fileschema.Variables(i).Dimensions(1).Length=L+1;
                fileschema.Variables(i).Dimensions(2).Name='eta_rho';
                fileschema.Variables(i).Dimensions(2).Length=M+1;
                fileschema.Variables(i).Dimensions(3).Name='tide_period';
                fileschema.Variables(i).Dimensions(3).Length=Nharmonics;
            end
            fileschema.Variables(i).Datatype='double';
            fileschema.Variables(i).Attributes(1).Name='long_name';
            fileschema.Variables(i).Attributes(1).Value=attvals{i}{1};
            fileschema.Variables(i).Attributes(2).Name='units';
            fileschema.Variables(i).Attributes(2).Value=attvals{i}{2};
        end
        disp('Creating netcdf file')
        ncwriteschema(fnOut,fileschema);
    end

%Write to output file
disp('Writing to netcdf file')
ncwrite(fnOut,'tide_period',ROMSperiods)
ncwrite(fnOut,'tide_Eamp',zamp)
ncwrite(fnOut,'tide_Ephase',zpha)
ncwrite(fnOut,'tide_Cmax',major)
ncwrite(fnOut,'tide_Cmin',minor)
ncwrite(fnOut,'tide_Cangle',inclination)
ncwrite(fnOut,'tide_Cphase',phase)
ncwrite(fnOut,'tide_Uamp',uamp)
ncwrite(fnOut,'tide_Uphase',upha)
ncwrite(fnOut,'tide_Vamp',vamp)
ncwrite(fnOut,'tide_Vphase',vpha)
end
%%
function [fFac,uFac]=TPXOnodalfactors(dnum,Names)

%f and u factors for the harmonics listed in cell array Names.
%(e.g. Names={'MM' 'M2' S2'};
%Only those which are in the TPXO model are evaluated.
%They are evaluated at time dnum (a Matlab datenumber).
%See Table xxvi in A.T. Doodson (1928) 'On the Analysis of Tidal Observations'
%Philosophical Transactions of the Royal Society of London. Series A, Vol. 227
%J. Luick, Thanksgiving Day, 2011, Adelaide

%if f and u are not reassigned below, they are probably solar
%terms, i.e. have f=1 and u=0.
f=struct;
u=struct;
% fFac=ones(length(Names),1);
% uFac=zeros(length(Names),1);

t=(dnum+0.5-datenum(1900,1,1))/36525;
VN=mod(360*(0.719954-5.372617*t+0.000006*t*t),360);
VN(VN<0)=VN(VN<0)+360;
VN=VN*pi/180;

%coefficients
cN=cos(VN);
c2N=cos(2*VN);
c3N=cos(3*VN);
sN=sin(VN);
s2N=sin(2*VN);
s3N=sin(3*VN);

%Assign values for f and u of nonsolar constituents
%Doodson Table XXVI (with u*pi/180)
f.MM=1.0-0.1300*cN+0.0013*c2N;
u.MM=0;
f.MF=1.0429+0.4135*cN-0.004*c2N;
u.MF=-0.4143*sN+0.0468*s2N-0.0066*s3N;

f.O1=1.0089+.1871*cN-0.0147*c2N+0.0014*c3N;
u.O1=0.1885*sN-0.0234*s2N+0.0033*s3N;

f.K1=1.0060+0.1150*cN-0.0088*c2N+0.0006*c3N;
u.K1=-0.1546*sN+0.0119*s2N-0.0012*s3N;

f.M2=1.0004-0.0373*cN+0.0002*c2N;
u.M2=-0.0374*sN;

f.K2=1.0241+0.2863*cN+0.0083*c2N-0.0015*c3N;
u.K2=-0.3096*sN+0.0119*s2N-0.0007*s3N;

% dependent values
f.Q1=f.O1;
u.Q1=u.O1;
f.N2=f.M2;
u.N2=u.M2;
f.MN4=f.M2^2;
u.MN4=2*u.M2;
f.M4=f.M2^2;
u.M4=2*u.M2;
f.MS4=f.M2;
u.MS4=u.M2;


%Assign fFac and uFac
for n=1:length(Names)
    if isfield(f,Names{n});
        fFac.(Names{n})=f.(Names{n});
        uFac.(Names{n})=mod(u.(Names{n})*180/pi,360);
    else
        %if f and u are not assigned, they are probably solar
        %terms, i.e. have f=1 and u=0.
        fFac.(Names{n})=1;
        uFac.(Names{n})=0;
    end
end

end


%%
% Compute equilibrium phases in accordance with Cartwright "tidal analysis - a
% retrospect", 1982, pp170 - 188 in "Time series methods in hydrosciences,
% A.H.el-Shaarawi and S.R.Esterby (eds), Elsevier
% dnum (Matlab datenumber) need not be an integer (V will be computed for
% the actual time, not the integral part).
%J. Luick, www.austides.com
%************************************************************
function Vdeg=Vphase(dnum,DN_List)
%************************************************************

t=(dnum+0.5-datenum(1900,1,1))/36525;
tHour=mod(dnum,1)*24;     %Hour of day

DN_Nbr=DN_List(1:7);
DN_Pha=DN_List(:,7);
DN_Spd=DN_List(:,8);

Vs=mod(360*(0.751206 + 1336.855231*t - 0.000003*t*t),360);
Vh=mod(360*(0.776935 +  100.002136*t + 0.000001*t*t),360);
Vp=mod(360*(0.928693 +   11.302872*t - 0.000029*t*t),360);
VN=mod(360*(0.719954 -    5.372617*t + 0.000006*t*t),360);
Vp1=mod(360*(0.781169 +   0.004775*t + 0.000001*t*t),360);
Vs(Vs<0)=Vs(Vs<0)+360;
Vh(Vs<0)=Vh(Vh<0)+360;
Vp(Vp<0)=Vp(Vp<0)+360;
VN(VN<0)=VN(VN<0)+360;
Vp1(Vp1<0)=Vp1(Vp1<0)+360;

Vdeg=tHour*DN_Spd+Vs*DN_Nbr(:,2)+Vh*DN_Nbr(:,3)+...
    Vp*DN_Nbr(:,4)+VN*DN_Nbr(:,5)+Vp1*DN_Nbr(:,6)+DN_Pha;
Vdeg=mod(Vdeg,360);
end
%%
function VARinterp=interpTPXO(TPXOvar,harmonic,count,lon,lat,mask)
%*****************************************************************
% Extract TPXO data and interpolate onto ROMS grid
% Inputs:
% VAR: hRE, hIm, uRE, uIm, vRe, or vIm
% harmonic: one of: M2 S2 N2 K2 K1 O1 P1 Q1 MF MM M4 MS4 MN4
% lon, lat: 2D arrays of longitude and latitude to interpolate to
% mask: array of 1s and 0s corresponding to wet and dry points of lon & lat
% Output: VARinterp (VAR on lon, lat grid)
% TPXO files must be on the matlab path
% J. Luick, www.austides.com
% Modified by C James October 2013
%*****************************************************************

VARinterp=zeros(size(mask));
iswet=mask==1;

if strcmp(harmonic,'MS4') || strcmp(harmonic,'MN4') || ...
      strcmp(harmonic,'MF') || strcmp(harmonic,'MM')
   x=TPXOvar.x_a6;
   y=TPXOvar.y_a6;
   z=squeeze(TPXOvar.z_a6(:,:,count));
   depth=TPXOvar.depth_a6;
   m=TPXOvar.mask_a6&depth>0;
else
   x=TPXOvar.x_a30;
   y=TPXOvar.y_a30;
   z=squeeze(TPXOvar.z_a30(:,:,count));
   depth=TPXOvar.depth_a30;
   m=TPXOvar.mask_a30&depth>0;
end

if(exist('scatteredInterpolant','file') && exist('griddedInterpolant','file'))
   F=scatteredInterpolant(x(m),y(m),z(m),'nearest');
   z=F(x,y);
   F=griddedInterpolant(x',y',z');
   V1=F(lon(iswet),lat(iswet))';    %z=F(x,y) interpolates to 1-d vectors x and y  (z also is 1-d)
   VARinterp(iswet)=V1;
else
   F=TriScatteredInterp(x(m),y(m),z(m),'nearest');
   z=F(x,y);
   F=TriScatteredInterp(x(:),y(:),z(:));
   V1=F(lon(iswet),lat(iswet));
   VARinterp(iswet)=V1;
end

end
%%
function [SEMA,  ECC, INC, PHA, w]=ap2ep(Au, PHIu, Av, PHIv, plot_demo)
% Convert tidal amplitude and phase lag (ap-) parameters into tidal ellipse
% (e-) parameters. Please refer to ep2app for its inverse function.
%
% Usage:
%
% [SEMA,  ECC, INC, PHA, w]=app2ep(Au, PHIu, Av, PHIv, plot_demo)
%
% where:
%
%     Au, PHIu, Av, PHIv are the amplitudes and phase lags (in degrees) of
%     u- and v- tidal current components. They can be vectors or
%     matrices or multidimensional arrays.
%
%     plot_demo is an optional argument, when it is supplied as an array
%     of indices, say [i j k l], the program will plot an  ellipse
%     corresponding to Au(i,j, k, l), PHIu(i,j,k,l), Av(i,j,k,l), and
%     PHIv(i,j,k,l);
%
%     Any number of dimensions are allowed as long as your computer
%     resource can handle.
%
%     SEMA: Semi-major axes, or the maximum speed;
%     ECC:  Eccentricity, the ratio of semi-minor axis over
%           the semi-major axis; its negative value indicates that the ellipse
%           is traversed in clockwise direction.
%     INC:  Inclination, the angles (in degrees) between the semi-major
%           axes and u-axis.
%     PHA:  Phase angles, the time (in angles and in degrees) when the
%           tidal currents reach their maximum speeds,  (i.e.
%           PHA=omega*tmax).
%
%           These four e-parameters will have the same dimensionality
%           (i.e., vectors, or matrices) as the input ap-parameters.
%
%     w:    Optional. If it is requested, it will be output as matrices
%           whose rows allow for plotting ellipses and whose columns are
%           for different ellipses corresponding columnwise to SEMA. For
%           example, plot(real(w(1,:)), imag(w(1,:))) will let you see
%           the first ellipse. You may need to use squeeze function when
%           w is a more than two dimensional array. See example.m.
%
% Document:   tidal_ellipse.ps

if nargin < 5
    plot_demo=0;  % by default, no plot for the ellipse
end


% Assume the input phase lags are in degrees and convert them in radians.
PHIu = PHIu/180*pi;
PHIv = PHIv/180*pi;

% Make complex amplitudes for u and v
i = sqrt(-1);
u = Au.*exp(-i*PHIu);
v = Av.*exp(-i*PHIv);

% Calculate complex radius of anticlockwise and clockwise circles:
wp = (u+i*v)/2;      % for anticlockwise circles
wm = conj(u-i*v)/2;  % for clockwise circles
% and their amplitudes and angles
Wp = abs(wp);
Wm = abs(wm);
THETAp = angle(wp);
THETAm = angle(wm);

% calculate e-parameters (ellipse parameters)
SEMA = Wp+Wm;              % Semi  Major Axis, or maximum speed
SEMI = Wp-Wm;              % Semin Minor Axis, or minimum speed
ECC = SEMI./SEMA;          % Eccentricity

PHA = (THETAm-THETAp)/2;   % Phase angle, the time (in angle) when
% the velocity reaches the maximum
INC = (THETAm+THETAp)/2;   % Inclination, the angle between the
% semi major axis and x-axis (or u-axis).

% convert to degrees for output
PHA = PHA/pi*180;
INC = INC/pi*180;
%THETAp = THETAp/pi*180;
%THETAm = THETAm/pi*180;

% flip THETAp and THETAm, PHA, and INC in the range of
% [-pi, 0) to [pi, 2*pi), which at least is my convention.
% id = THETAp < 0;   THETAp(id) = THETAp(id)+360;
% id = THETAm < 0;   THETAm(id) = THETAm(id)+360;
id = PHA < 0;      PHA(id) = PHA(id)+360;
id = INC < 0;      INC(id) = INC(id)+360;


if nargout == 5
    ndot=36;
    dot=2*pi/ndot;
    ot=0:dot:2*pi-dot;
    w=wp(:)*exp(i*ot)+wm(:)*exp(-i*ot);
    w=reshape(w, [size(Au) ndot]);
end


if any(plot_demo)
    plot_ell(SEMA, ECC, INC, PHA, plot_demo)
end
end

%%
%Authorship Copyright:
%
%    The author of this program retains the copyright of this program, while
% you are welcome to use and distribute this program as long as you credit
% the author properly and respect the program name itself. Particularly,
% you are expected to retain the original author's name in this original
% version of the program or any of its modified version that you might make.
% You are also expected not to essentially change the name of the programs
% except for adding possible extension for your own version you might create,
% e.g. app2ep_xx is acceptable.  Any suggestions are welcome and enjoying my
% program(s)!
%
%
%Author Info:
%_______________________________________________________________________
%  Zhigang Xu, Ph.D.
%  (pronounced as Tsi Gahng Hsu)
%  Research Scientist
%  Coastal Circulation
%  Bedford Institute of Oceanography
%  1 Challenge Dr.
%  P.O. Box 1006                    Phone  (902) 426-2307 (o)
%  Dartmouth, Nova Scotia           Fax    (902) 426-7827
%  CANADA B2Y 4A2                   email zhigangx@emerald.bio.dfo.ca
%                                         zhigang_xu_98@yahoo.com
%_______________________________________________________________________
%
%Release Date: Nov. 2000
%%
function TPXO=readTPXOdata(TPXOfile,ROMSnames,lon,lat)
%Read TPXO version 8 data
%J. Luick, www.austides.com
%With help from C. James

lonmin=min(min(lon));
lonmax=max(max(lon));
latmin=min(min(lat));
latmax=max(max(lat));

bndx=[lon(1,1) lon(end,1) lon(end,end) lon(1,end)];
bndy=[lat(1,1) lat(end,1) lat(end,end) lat(1,end)];

var={'h','U','V'};

count_a6=0;
count_a30=0;
for n=1:length(ROMSnames)
   str=['Harmonic: ' ROMSnames(n)];
   disp(str)
   if strcmp(ROMSnames{n},'MS4') || strcmp(ROMSnames{n},'MN4') || ...
         strcmp(ROMSnames{n},'MF') || strcmp(ROMSnames{n},'MM')
      TPXOfile_grid=TPXOfile.grid.atlas6;
      count_a6=count_a6+1;
      TPXO.('harmonic_a6').harmonic(count_a6)=ROMSnames(n);
   else
      TPXOfile_grid=TPXOfile.grid.atlas30;
      count_a30=count_a30+1;
      TPXO.('harmonic_a30').harmonic(count_a30)=ROMSnames(n);
   end
   TPXOfile_elev=TPXOfile.elev.(lower(ROMSnames{n}));
   TPXOfile_vel=TPXOfile.vel.(lower(ROMSnames{n}));
   file={TPXOfile_elev,TPXOfile_vel,TPXOfile_vel};
   for i=1:length(var);
      if strcmp(var{i},'h')
         coord='z';
      elseif strcmp(var{i},'U')
         coord='u';
      elseif strcmp(var{i},'V')
         coord='v';
      end
      X=ncread(file{i},['lon_' coord]);
      Y=ncread(file{i},['lat_' coord]);
      I=find((Y>=latmin-0.5)&(Y<=latmax+0.5));
      J=find((X>=lonmin-0.5)&(X<=lonmax+0.5));
      istart=I(1);
      jstart=J(1);
      icount=length(I);
      jcount=length(J);

      H=ncread(TPXOfile_grid,['h' coord],[istart jstart],[icount jcount]);
      Z=complex(ncread(file{i},[lower(var{i}) 'Re'],[istart jstart],[icount jcount]),...
                ncread(file{i},[lower(var{i}) 'Im'],[istart jstart],[icount jcount]));

      [x,y]=meshgrid(X(J),Y(I));
      if strcmp(var{i},'h')
         z=double(Z)/1000;
      elseif strcmp(var{i},'U') || strcmp(var{i},'V')
         z=double(Z)/10000;
      end
      depth=double(H);
      mask=inpolygon(x,y,bndx,bndy);
      if strcmp(TPXOfile_grid,'grid_tpxo8_atlas6.nc');
         TPXO.(var{i}).x_a6=x;
         TPXO.(var{i}).y_a6=y;
         TPXO.(var{i}).depth_a6=depth;
         TPXO.(var{i}).mask_a6=mask;
         if strcmp(var{i},'h');
            TPXO.(var{i}).z_a6(:,:,count_a6)=z;
         else
            TPXO.(var{i}).z_a6(:,:,count_a6)=z./repmat(TPXO.(var{i}).depth_a6,[1 1 size(z,3)]);
         end
      elseif strcmp(TPXOfile_grid,'grid_tpxo8atlas_30.nc');
         TPXO.(var{i}).x_a30=x;
         TPXO.(var{i}).y_a30=y;
         TPXO.(var{i}).depth_a30=depth;
         TPXO.(var{i}).mask_a30=mask;
         if strcmp(var{i},'h');
            TPXO.(var{i}).z_a30(:,:,count_a30)=z;
         else
            TPXO.(var{i}).z_a30(:,:,count_a30)=z./repmat(TPXO.(var{i}).depth_a30,[1 1 size(z,3)]);
         end
      end
   end
end

end

%%
function varcheck(varname,default_value)
% function varcheck(varname,default_value)
% checks calling workspace for existence of variable var
% if it is non-existent it creates it in gives it the value default_value
% if it exists but is empty it also assigns it the value default_value
% if it exists and has any other value it is not altered.
% useful for testing optional inputs into functions
% Charles James 2012

if (nargin<2)||~ischar(varname);
    return;
end

a=evalin('caller',['exist(''' varname ''',''var'');']);

if (a~=0)
    var=evalin('caller',varname);
    if isempty(var)
        var=default_value;
    end
else
    var=default_value;
end

assignin('caller',varname,var);

end
