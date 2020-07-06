function aaa___test___ETA_disp(varargin)
%no actual test is required, because there is no input checking

%if nargin==0,RunTestHeadless=false;else,RunTestHeadless=true;end

% suppress v6.5 warning
v=version;v(strfind(v,'.'):end)='';v=str2double(v);
if v<=7
    warning off MATLAB:m_warning_end_without_block;clc
end

try
    N=10;then=now;
    for n=1:N
        pause((1+rand)*2/N)
        ETA_disp(N,n,then)
    end
    % N=10^6;then=now;
    % for n=1:N
    %     pause((1+rand)*2/N)
    %     if mod(n,10^4)==0 %avoid spending too much time printing the ETA
    %         ETA_disp(N,n,then)
    %     end
    % end
catch
    error('examples failed')
end

then=now-0.001;
s=ETA_disp(2,1,then);
s(2)=[];%depends on current time, so remove it before hashing
%ComputeNonCryptHash(s,32),return
if ~strcmp(ComputeNonCryptHash(s,32),'DB8302F6')
    error('hash didn''t match (so either the hash function changed, or the ETA_disp output)')
else
    disp('test completed succesfully')
end
end
function data=cast_to_uint16_vector(data)
%Linearize the input data and convert it to a uint16 vector
data=cast_to_uint16_vector__cell({data});
data([end-1 end])=[];%Remove the [1 1] that is added because of the wrapping in a cell
end
function data=cast_to_uint16_vector__cell(data)
sz=size(data);data=data(:);
for n=1:numel(data)
    if numel(data{n})==0
        c=double(class(data{n})');
        data{n}=uint16([0;c;size(data{n})']);
        continue
    end
    switch class(data{n})
        case {'double','single'}
            data{n}=cast_to_uint16_vector__floats(data{n});
        case 'logical'
            data{n}=cast_to_uint16_vector__logical(data{n});
        case {'uint8','uint16','uint32','uint64','int8','int16','int32','int64'}
            data{n}=cast_to_uint16_vector__integer(data{n});
        case 'char'
            data{n}=cast_to_uint16_vector__char(data{n});
        case 'cell'
            data{n}=cast_to_uint16_vector__cell(data{n});
        case 'struct'
            data{n}=cast_to_uint16_vector__struct(data{n});
        otherwise
            error('HJW:cast_to_uint16_vector:nosupport',...
                'Unsupported data type in nested variable')
    end
end
data=cell2mat(data);%Merge all cell contents
data=[data;uint16(mod(sz',2^16))];%Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__floats(data)
sz=size(data);c=class(data);%the rest of the function treats singles as double

%Convert to a uint64, separate it into 4 words and everything merge into a vector
[bit64,bit4]=typecast_double_uint64(double(data));
bit4_round=mod(bit64,2^16);bit64=bit64-bit4_round;bit64=bit64/2^16;
bit3      =mod(bit64,2^16);bit64=bit64-bit3;      bit64=bit64/2^16;
bit2      =mod(bit64,2^16);bit64=bit64-bit2;      bit64=bit64/2^16;
bit1      =mod(bit64,2^16);
data=[bit1';bit2';bit3';bit4'];
data=uint16(data(:));

%Append the array size and type to make it influence the hash.
data=[data;uint16(c.');uint16(mod(sz',2^16))];
end
function data=cast_to_uint16_vector__logical(data)
sz=size(data);data=data(:);

if mod(numel(data),16) %pad to 16 bits with 0
    data(16*ceil(numel(data)/16))=0;
end
vector=uint16(2.^(15:-1:0))';
data=uint16(reshape(data,16,[]));
try
    data=data.*vector;
catch %no implicit expansion
    %times() is not defined for uint16 on ML6.5
    data=double(data).*repmat(double(vector),[1 size(data,2)]);
    data=uint16(data);
end
data=uint16(sum(data,1)).';

data=[data;uint16(mod(sz',2^16))];%Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__integer(data)
%Large values (>2^52) will not have integer precision due to a conversion to double.
%This conversion is done, because of limited availability of operations on ML6.5.
sz=size(data);data=data(:);

c=class(data);
if c(1)~='u'
    %Shift int* values to the uint* data range
    data=double(data)-double(eval([c '(-inf)']));
else
    data=double(data);
end
switch c(end)
    case '8'
        %append a 0 for odd length and merge pairs
        if mod(numel(data),2),data(end+1)=0;end
        data=reshape(data,[],2);
        data=data(:,1)*255+data(:,2);
        data=uint16(data);
    case '6'
        data=uint16(data);
    case '2'
        %split to 2 words
        bit1=floor(data/2^16);
        bit2=mod(data,2^16);
        data=[bit1';bit2'];
        data=uint16(data(:));
    case '4'
        %split to 4 words, bit4 contains a rounding error for data >2^52
        bit64=data;
        bit4_round=mod(bit64,2^16);bit64=bit64-bit4_round;bit64=bit64/2^16;
        bit3      =mod(bit64,2^16);bit64=bit64-bit3;      bit64=bit64/2^16;
        bit2      =mod(bit64,2^16);bit64=bit64-bit2;      bit64=bit64/2^16;
        bit1      =mod(bit64,2^16);
        
        data=[bit1';bit2';bit3';bit4_round'];
        data=uint16(data(:));
end
%Append the array size and type to make it influence the hash.
data=[data;uint16(c.');uint16(mod(sz',2^16))];
end
function data=cast_to_uint16_vector__char(data)
sz=size(data);data=data(:);
data=uint16(data);%Normal chars are already a uint16 internally (in Matlab).
data=[data;uint16(mod(sz',2^16))];%Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__struct(data)
sz=size(data);data=data(:);
fn=fieldnames(data);
output=cell(2,numel(fn));
for n=1:numel(fn)
    output{1,n}=fn{n};
    output{2,n}={data.(fn{n})};
end
data=cast_to_uint16_vector__cell(output);
data=[data;uint16(mod(sz',2^16))];%Append the array size to make it influence the hash.
end
function hash=ComputeNonCryptHash(data,varargin)
%Compute a non-cryptographic hash
%
% This function is intended to be fast, but without requiring a Java or mex implementation to do
% the actual hashing. It was *not* checked for any security flaws and is therefore probably
% vulnerable to most attacks.
% Non-cryptographic hashes should only be used as a checksum. Don't use this to do things like
% storing passwords.
%
%syntax:
%  hash=ComputeNonCryptHash(data)
%  hash=ComputeNonCryptHash(data,HashLength)
%
%data       The data to be hashed. Most common data types are allowed: uint*, int*, char, cell,
%           struct, double, or single (string is cast to char). The contents of the nested data
%           types (i.e. cell and struct) must also be one of the mentioned data types.
%HashLength The length of the hash (the number of bits). This value must be a multiple of 16. The
%           default is 256 bits. Depending on your input 64 bits might have some collisions, but
%           64 bits and higher should be safe.
%
%hash       The hash in an upper case hexadecimal char vector of size 1x(HashLength/4).
%
%  _______________________________________________________________________
% | Compatibility | Windows 10  | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |---------------|-------------|------------------|----------------------|
% | ML R2020a     |  works      |  not tested      |  not tested          |
% | ML R2018a     |  works      |  works           |  not tested          |
% | ML R2015a     |  works      |  works           |  not tested          |
% | ML R2011a     |  works      |  works           |  not tested          |
% | ML 6.5 (R13)  |  works      |  not tested      |  not tested          |
% | Octave 5.2.0  |  works      |  works           |  not tested          |
% | Octave 4.4.1  |  works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
%
% Version: 1.0
% Date:    2020-06-11
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( http://creativecommons.org/licenses/by-nc-sa/4.0 )
% Email=  'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

if nargin<1
    error('HJW:ComputeNonCryptHash:InputIncorrect','At least 1 input required.')
else
    try
        if isa(data,'string'),data=char(data);end
    catch
        error('HJW:ComputeNonCryptHash:ConvertToChar',...
            'The required conversion from string to char failed.')
    end
end
if nargin<2
    HashLength=256;
else
    HashLength=varargin{1};
    if numel(HashLength)~=1 || ~isnumeric(HashLength) || mod(HashLength,16)~=0 || HashLength<16
        error('HJW:ComputeNonCryptHash:InputIncorrect',...
            'Second input (hash length) must be a multiple of 16.')
    end
end

try
    %Convert the input to an uint16 array (Nx1).
    data=cast_to_uint16_vector(data);
catch
    ME=lasterror; %#ok<LERR>
    if strcmp(ME.identifier,'MATLAB:nomem')
        %rethrow memory error
        rethrow(ME)
    else
        error('HJW:ComputeNonCryptHash:UnwindFailed',...
            'The nested input contains an unsupported data type.')
    end
end

%Extend to a multiple of HashLength bits. Padding with zeros is generally not advised, and the
%performance penalty for this extension (compared to padding with zeros) should be negligible.
if mod(numel(data),HashLength/16)
    extra=uint16(1:HashLength/16).'; extra(1:mod(numel(data),HashLength/16))=[];
    data=[data;extra];
end

%Add perturbation to the data and convert to 16xN logical. Then further perturb the intermediate
%result by doing an XOR in blocks with itself (by reshaping and transposing).
data=ComputeNonCryptHash_shuffle_uint16(data);
data=ComputeNonCryptHash_uint16_to_logical(data);
data=xor(data,reshape(data,[],16).');

%Reshape to HashLength cols and collapse the key size down to the hash length by counting the
%number of true bits (even=1, odd=0).
data=mod(sum(reshape(data,HashLength,[]),2),2);
data=ComputeNonCryptHash_logical_to_uint16(data);

if nargin>2
    hash=data;%Return uint16 for the salting.
    return
end

%Perturb the hash, analogous to salting. This function computes the hash of the hash and applies a
%few operations to the data to increase the randomness of the end result.
data=ComputeNonCryptHash_add_salt(data);

%Convert the (HashLength/16)x1 uint16 to a hash string by encoding it as hexadecimal.
hash=ComputeNonCryptHash_dec2hex(data);hash=reshape(hash.',1,[]);
end
function data=ComputeNonCryptHash_add_salt(data)
%Apply a few transformations to the hash to increase the spread.
%If this function is not added, the hashes of -12345678 and 12345678 will be very similar.
%A true salt would be to append the salt to the data, so this is not actually a salt.
saltHashLength=16*numel(data);
salt=ComputeNonCryptHash(data,saltHashLength,[]);%avoid infinite recursion by using a third input
salt=ComputeNonCryptHash_shuffle_uint16_inv(salt);
data=mod(double(data).*double(salt),1+2^16);
data=uint16(data);
end
function hash=ComputeNonCryptHash_dec2hex(data)
%Look up the precomputed dec2hex for faster conversion.
persistent LUT
if isempty(LUT)
    LUT=upper(dec2hex(0:(-1+2^16),4));%even though the default should already upper case
end
data=double(data)+1;% plus() is not defined for uint16 on ML6.5
hash=LUT(data,:);
end
function data=ComputeNonCryptHash_logical_to_uint16(data)
if mod(numel(data),16) %pad to 16 bits with 0
    data(16*ceil(numel(data)/16))=0;
end
vector=uint16(2.^(15:-1:0))';
data=uint16(reshape(data,16,[]));
try
    data=data.*vector;
catch %no implicit expansion
    %times() is not defined for uint16 on ML6.5
    data=double(data).*repmat(double(vector),[1 size(data,2)]);
    data=uint16(data);
end
data=uint16(sum(data,1)).';
end
function data=ComputeNonCryptHash_shuffle_uint16(data)
%input should be uint16
base=65537;%base=(1+(2^16));
key=479001600;%key=1*2*3*4*5*6*7*8*9*10*11*12;
data = uint16(mod(double(data) * key , base));
end
function data=ComputeNonCryptHash_shuffle_uint16_inv(data)
base=65537;%base=(1+(2^16));
%key=1*2*3*4*5*6*7*8*9*10*11*12;
% %solution suggested by John D'Errico, https://www.mathworks.com/matlabcentral/answers/81859
% [G,C]=gcd(key,base);invKey=mod(C,base);
invKey=1919;
data=uint16(mod(double(data) * invKey,base));
end
function data=ComputeNonCryptHash_uint16_to_logical(data)
%uint16 Nx1 vector in, logical 16xN array out
persistent LUT
if isempty(LUT)
    LUT=dec2bin(0:(-1+2^16))=='1';
    LUT=LUT.';
end
data=double(data)+1;% plus() is not defined for uint16 on ML6.5
data=LUT(:,data);
end
function out=ETA_disp(N_total,n,then)
%Display the estimated time of completion of a looping process.
%
%input:
%  N_total: total number of iterations
%  n: current count of iterations
%  then: the output of the now function when the loop began
%output:
%  if called without output arguments, the message is printed to the command window, otherwise the
%  message is returned in a cell array.
%
% The output for tic was introduced in R2008b, so tic/toc can't be nested in all releases, while
% this implementation can be nested in any release.
%
% To keep this function fast no input checking is done.
%
% Examples:
%
%   N=10;then=now;
%   for n=1:N
%       pause((1+rand)*2/N)
%       ETA_disp(N,n,then)
%   end
%   
%   N=10^6;then=now;
%   for n=1:N
%       pause((1+rand)*2/N)
%       if mod(n,10^4)==0 %avoid spending too much time printing the ETA
%           ETA_disp(N,n,then)
%       end
%   end
%
%  _______________________________________________________________________
% | Compatibility | Windows 10  | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |---------------|-------------|------------------|----------------------|
% | ML R2020a     |  works      |  not tested      |  not tested          |
% | ML R2018a     |  works      |  works           |  not tested          |
% | ML R2015a     |  works      |  works           |  not tested          |
% | ML R2011a     |  works      |  works           |  not tested          |
% | ML 6.5 (R13)  |  works      |  not tested      |  not tested          |
% | Octave 5.2.0  |  works      |  works           |  not tested          |
% | Octave 4.4.1  |  works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
%
% Version: 1.0.1
% Date:    2020-06-11
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( http://creativecommons.org/licenses/by-nc-sa/4.0 )
% Email=  'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

clc
persistent legacy
if isempty(legacy)
    %addtodate was probably added in ML7.0
    legacy.addtodate=ifversion('<',7,'Octave','<',4);
end

t=(now-then)*(24*60*60);%equivalent to t=toc(h_tic);
t_min_total=floor(t/60);t_sec_total=round(t-60*t_min_total);
f=n/N_total;t_r=round((t-f*t)/f);
if legacy.addtodate
    ETA=datestr(now+(t_r/(24*60*60)),'HH:MM:SS');
else
    ETA=datestr(addtodate(now,t_r,'second'),'HH:MM:SS');
end
t_per_it=t/n;
t_min=floor(t_per_it/60);t_sec=t_per_it-60*t_min;

out=cell(1,3);
out{1}=sprintf('%05.1f%% done after a total time of %02d:%02d.\n',...
    100*n/N_total,t_min_total,t_sec_total);
out{2}=sprintf('estimated time of completion: %s\n',...
    ETA);
out{3}=sprintf('(%02d:%05.2f per iteration, %d iterations left)\n',...
    t_min,t_sec,N_total-n);
if nargout==0
    fprintf('%s',out{:});
    clear out;
end
end
function tf=ifversion(test,Rxxxxab,Oct_flag,Oct_test,Oct_ver)
%Determine if the current version satisfies a version restriction
%
% To keep the function fast, no input checking is done. This function returns a NaN if a release
% name is used that is not in the dictionary.
%
% Syntax:
% tf=ifversion(test,Rxxxxab)
% tf=ifversion(test,Rxxxxab,'Octave',test_for_Octave,v_Octave)
%
% Output:
% tf       - If the current version satisfies the test this returns true.
%            This works similar to verLessThan.
%
% Inputs:
% Rxxxxab - Char array containing a release description (e.g. 'R13', 'R14SP2' or 'R2019a') or the
%           numeric version.
% test    - Char array containing a logical test. The interpretation of this is equivalent to
%           eval([current test Rxxxxab]). For examples, see below.
%
% Examples:
% ifversion('>=','R2009a') returns true when run on R2009a or later
% ifversion('<','R2016a') returns true when run on R2015b or older
% ifversion('==','R2018a') returns true only when run on R2018a
% ifversion('==',9.8) returns true only when run on R2020a
% ifversion('<',0,'Octave','>',0) returns true only on Octave
%
% The conversion is based on a manual list and therefore needs to be updated manually, so it might
% not be complete. Although it should be possible to load the list from Wikipedia, this is not
% implemented.
%
%  _______________________________________________________________________
% | Compatibility | Windows 10  | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |---------------|-------------|------------------|----------------------|
% | ML R2020a     |  works      |  not tested      |  not tested          |
% | ML R2018a     |  works      |  works           |  not tested          |
% | ML R2015a     |  works      |  works           |  not tested          |
% | ML R2011a     |  works      |  works           |  not tested          |
% | ML 6.5 (R13)  |  works      |  not tested      |  not tested          |
% | Octave 5.2.0  |  works      |  works           |  not tested          |
% | Octave 4.4.1  |  works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
%
% Version: 1.0.2
% Date:    2020-05-20
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 )
% Email=  'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

%The decimal of the version numbers are padded with a 0 to make sure v7.10 is larger than v7.9.
%This does mean that any numeric version input needs to be adapted. multiply by 100 and round to
%remove the potential for float rounding errors.
%Store in persistent for fast recall (don't use getpref, as that is slower than generating the
%variables and makes updating this function harder).
persistent  v_num v_dict octave
if isempty(v_num)
    %test if Octave is used instead of Matlab
    octave=exist('OCTAVE_VERSION', 'builtin');
    
    %get current version number
    v_num=version;
    ii=strfind(v_num,'.');
    if numel(ii)~=1,v_num(ii(2):end)='';ii=ii(1);end
    v_num=[str2double(v_num(1:(ii-1))) str2double(v_num((ii+1):end))];
    v_num=v_num(1)+v_num(2)/100;
    v_num=round(100*v_num);%remove float rounding errors
    
    %get dictionary to use for ismember
    v_dict={...
        'R13' 605;'R13SP1' 605;'R13SP2' 605;'R14' 700;'R14SP1' 700;'R14SP2' 700;'R14SP3' 701;...
        'R2006a' 702;'R2006b' 703;'R2007a' 704;'R2007b' 705;'R2008a' 706;'R2008b' 707;...
        'R2009a' 708;'R2009b' 709;'R2010a' 710;'R2010b' 711;'R2011a' 712;'R2011b' 713;...
        'R2012a' 714;'R2012b' 800;'R2013a' 801;'R2013b' 802;'R2014a' 803;'R2014b' 804;...
        'R2015a' 805;'R2015b' 806;'R2016a' 900;'R2016b' 901;'R2017a' 902;'R2017b' 903;...
        'R2018a' 904;'R2018b' 905;'R2019a' 906;'R2019b' 907;'R2020a' 908};
end

if octave
    if nargin==2
        warning('HJW:ifversion:NoOctaveTest',...
            ['No version test for Octave was provided.',char(10),...
            'This function might return an unexpected outcome.']) %#ok<CHARTEN>
        %Use the same test as for Matlab, which will probably fail.
        L=ismember(v_dict(:,1),Rxxxxab);
        if sum(L)~=1
            warning('HJW:ifversion:NotInDict',...
                'The requested version is not in the hard-coded list.')
            tf=NaN;return
        else
            v=v_dict{L,2};
        end
    elseif nargin==4
        %undocumented shorthand syntax: skip the 'Octave' argument
        [test,v]=deal(Oct_flag,Oct_test);
        %convert 4.1 to 401
        v=0.1*v+0.9*fix(v);v=round(100*v);
    else
        [test,v]=deal(Oct_test,Oct_ver);
        %convert 4.1 to 401
        v=0.1*v+0.9*fix(v);v=round(100*v);
    end
else
    %convert R notation to numeric and convert 9.1 to 901
    if isnumeric(Rxxxxab)
        v=0.1*Rxxxxab+0.9*fix(Rxxxxab);v=round(100*v);
    else
        L=ismember(v_dict(:,1),Rxxxxab);
        if sum(L)~=1
            warning('HJW:ifversion:NotInDict',...
                'The requested version is not in the hard-coded list.')
            tf=NaN;return
        else
            v=v_dict{L,2};
        end
    end
end
switch test
    case '=='
        tf= v_num == v;
    case '<'
        tf= v_num <  v;
    case '<='
        tf= v_num <= v;
    case '>'
        tf= v_num >  v;
    case '>='
        tf= v_num >= v;
end
end
function [bit64,last16bits]=typecast_double_uint64(FP)
%Turn a double into a uint64 with the same binary representation.
%This is similar to typecast(FP,'uint64'); the difference being that this function ignores
%endianness, supports array inputs, and is slower.
%
%Because of missing support for some operations in ML6.5, this function returns the uint64 as a
%double. Because this may cause rounding errors, the last 16 bits are returned separately.

[M,E]=log2(FP);
signBit =-floor(sign(FP)/2-0.5);
exponent=E+1022;
mantissa=abs(M)*2-1;

%no plus() for integer types in ML6.5, so we need to use double, instead of uint64
bit64=zeros(size(FP));
bit64=bit64+(signBit*2^63);
bit64=bit64+(exponent*2^52);
bit64=bit64+(mantissa*2^52);
last16bits=mod(mantissa*2^52,2^16);

%correct the 0 and hard-code the special cases
L=isinf(FP);
bit64(FP==0)=0;
bit64(isnan(FP))=18444492273895866368;
bit64(L & FP>0)=9218868437227405312;%positive inf
bit64(L & FP<0)=18442240474082181120;%negative inf
last16bits(FP==0)=0;
last16bits(isnan(FP))=0;
last16bits(L)=0;
end