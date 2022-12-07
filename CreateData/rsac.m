%RSAC    Read SAC binary files.
%    RSAC('sacfile') reads in a SAC (seismic analysis code) binary
%    format file into a 3-column vector.
%    Column 1 contains time values.
%    Column 2 contains amplitude values.
%    Column 3 contains all SAC header information.
%    Default byte order is big-endian.  M-file can be set to default
%    little-endian byte order.
%
%    usage:  output = rsac('sacfile')
%
%    Examples:
%
%    KATH = rsac('KATH.R');
%    plot(KATH(:,1),KATH(:,2))
%
%    [SQRL, AAK] = rsac('SQRL.R','AAK.R');
%
%    by Michael Thorne (4/2004)   mthorne@asu.edu

function [varargout] = rsac(varargin);

for nrecs = 1:nargin

  sacfile = varargin{nrecs};

%---------------------------------------------------------------------------
%    Default byte-order
%    endian  = 'big'  big-endian byte order (e.g., UNIX)
%            = 'lil'  little-endian byte order (e.g., LINUX)

endian = 'lil';

if endian == 'big'
  fid = fopen(sacfile,'r','ieee-be'); 
elseif endian == 'lil'
  fid = fopen(sacfile,'r','ieee-le'); 
end

% read in single precision real header variables:
%---------------------------------------------------------------------------
for i=1:70
  h(i) = fread(fid,1,'single');
end


% read in single precision integer header variables:
%---------------------------------------------------------------------------
for i=71:105
  h(i) = fread(fid,1,'int32');
end

% read in logical header variables
%---------------------------------------------------------------------------
for i=106:110
  h(i) = fread(fid,1,'int32');
end

% read in character header variables
%---------------------------------------------------------------------------
for i=111:302
  h(i) = (fread(fid,1,'char'))';
end

% read in amplitudes
%---------------------------------------------------------------------------

YARRAY = fread(fid,'single');
ylen=length(YARRAY);

XARRAY = (linspace(h(6),(h(6)+ylen-1)*h(1),ylen))'; 



% add header signature for testing files for SAC format
%---------------------------------------------------------------------------
h(303) = 77;
h(304) = 73;
h(305) = 75;
h(306) = 69;

% arrange output files
%---------------------------------------------------------------------------


OUTPUT(:,1) = XARRAY;
OUTPUT(:,2) = YARRAY;

OUTPUT(1:306,3) = h(1:306)';

%pad xarray and yarray with NaN if smaller than header field
if h(80) < 306
  OUTPUT((h(80)+1):306,1) = NaN;
  OUTPUT((h(80)+1):306,2) = NaN;
end

fclose(fid);

varargout{nrecs} = OUTPUT;

end
