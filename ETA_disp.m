function [out,out_struct]=ETA_disp(N_total,n,then,per_it_display,scale_estimate)
%Display the estimated time of completion of a looping process
%
% To keep this function fast no input checking is done.
% The time remaining is capped at 10^10 seconds.
% The output for tic was introduced in R2008b, so tic/toc can't be nested in all releases, while
% this implementation can be nested in any release.
%
% Syntax:
%   ETA_disp(N_total,n,then)
%   ETA_disp(N_total,n,then,per_it_display)
%   ETA_disp(N_total,n,then,per_it_display,scale_estimate)
%   msg = ETA_disp(___)
%   [msg,strct] = ETA_disp(___)
%
% Input/output arguments:
% msg:
%   If an output argument is provided the message is not printed to the command window, but is
%   instead returned as a cell array.
% strct:
%   The second output argument is a struct with several fields containing the components of the
%   message. The included variables are the percentage done ('p_done'), the text with the elapsed
%   time ('elapsed_str'), the estimated time of completion ('ETA'), the estimated date of
%   completion ('ETA_date', empty if the estimated time left is less than a day), the number of
%   minutes per iteration ('t_min'), the number of seconds per iteration ('t_sec'), and the number
%   of iterations per second ('it_per_sec').
% N_total:
%   The total number of iterations.
% n:
%   The current count of iterations (this will be adjusted up to eps if 0).
% then:
%   This should be the value the now function returned before the loop began.
% per_it_display:
%   (optional, default is 0) 0 for automatic, -1 to display time per iteration on the last line, or
%   1 to print iterations per seconds on the last line.
% scale_estimate:
%   the time remaining is multiplied by this factor (this doesn't change the rate calculated over
%   the processed iterations). This allows skewing the estimate in cases where the earlier
%   iterations are expected to be faster.
%
% Examples:
%
%   N = 10;then = now; %#ok<TNOW1>
%   for n=1:N
%       pause((1+rand)*2/N)
%       ETA_disp(N,n,then)
%   end
%
%   N = 10^6;then = now; %#ok<TNOW1>
%   for n=1:N
%       pause((1+rand)*2/N)
%       if mod(n,10^4)==0 % Avoid spending too much time printing the ETA.
%           ETA_disp(N,n,then,1) % Print iterations/sec instead of sec/iteration.
%       end
%   end
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 2.0.0                                                         |%
%|  Date:    2023-11-03                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). You can see the full test matrix below.
% Compatibility considerations:
% - This is expected to work on all releases.
%
% /=========================================================================================\
% ||                     | Windows             | Linux               | MacOS               ||
% ||---------------------------------------------------------------------------------------||
% || Matlab R2023b       | W11: Pass           | Ubuntu 22.04: Pass  | Monterey: Pass      ||
% || Matlab R2023a       | W11: Pass           |                     |                     ||
% || Matlab R2022b       | W11: Pass           | Ubuntu 22.04: Pass  | Monterey: Pass      ||
% || Matlab R2022a       | W11: Pass           |                     |                     ||
% || Matlab R2021b       | W11: Pass           | Ubuntu 22.04: Pass  | Monterey: Pass      ||
% || Matlab R2021a       | W11: Pass           |                     |                     ||
% || Matlab R2020b       | W11: Pass           | Ubuntu 22.04: Pass  | Monterey: Pass      ||
% || Matlab R2020a       | W11: Pass           |                     |                     ||
% || Matlab R2019b       | W11: Pass           | Ubuntu 22.04: Pass  | Monterey: Pass      ||
% || Matlab R2019a       | W11: Pass           |                     |                     ||
% || Matlab R2018a       | W11: Pass           | Ubuntu 22.04: Pass  |                     ||
% || Matlab R2017b       | W11: Pass           | Ubuntu 22.04: Pass  | Monterey: Pass      ||
% || Matlab R2016b       | W11: Pass           | Ubuntu 22.04: Pass  | Monterey: Pass      ||
% || Matlab R2015a       | W11: Pass           | Ubuntu 22.04: Pass  |                     ||
% || Matlab R2013b       | W11: Pass           |                     |                     ||
% || Matlab R2007b       | W11: Pass           |                     |                     ||
% || Matlab 6.5 (R13)    | W11: Pass           |                     |                     ||
% || Octave 8.2.0        | W11: Pass           |                     |                     ||
% || Octave 7.2.0        | W11: Pass           |                     |                     ||
% || Octave 6.2.0        | W11: Pass           | Raspbian 11: Pass   | Catalina: Pass      ||
% || Octave 5.2.0        | W11: Pass           |                     |                     ||
% || Octave 4.4.1        | W11: Pass           |                     | Catalina: Pass      ||
% \=========================================================================================/

% Before doing any processing, make sure the time is captured to avoid any offsets caused by this
% function itself.
t = (now-then)*(24*60*60); %#ok<TNOW1> This is equivalent to t=toc(h_tic);

persistent legacy
if isempty(legacy)
    % The addtodate function was introduced in ML7.0 (R14), but the 'second' option was only
    % implemented later, in R2008b.
    legacy.addtodate = ifversion('<','R2008b','Octave','<',4);
end
if nargin<4 || numel(per_it_display)~=1,per_it_display = 0;end
if nargin<5 || numel(scale_estimate)~=1,scale_estimate = 1;end
n = max(eps,n); % Ensure n is non-zero, which would cause problems later.
t_min_total = floor(t/60);t_sec_total=round(t-60*t_min_total); % Calculate time elapsed.
f = n/N_total;t_r = (t-f*t)/f; % Calculation fraction done and time remaining.
if nargin>=5,t_r = t_r*scale_estimate;end % Scale the remainder estimate.
t_r = min(10^10,round(t_r)); % Round and cap to 10^10 seconds.
if legacy.addtodate
    ETA = datestr(now+(t_r/(24*60*60)),'HH:MM:SS'); %#ok<TNOW1,DATST>
else
    ETA = datestr(addtodate(now,t_r,'second'),'HH:MM:SS'); %#ok<DATST,DATOD,TNOW1>
end
days_left = floor(t_r/(24*60*60));
if days_left>0
    ETA_date = datestr(now+days_left,26); %#ok<TNOW1,DATST>
    ETA_ = sprintf('%s (on %s)',ETA,ETA_date);
else
    ETA_date = '';
    ETA_ = ETA;
end
t_per_it = t/n; % This is in seconds per iteration.
t_min = floor(t_per_it/60);t_sec=t_per_it-60*t_min;
if abs(per_it_display)~=1 % Revert to dynamic if the input is invalid.
    if t_per_it>0.5
        per_it_display = -1; % Display time per iteration.
    else
        per_it_display =  1; % Display iterations per second.
    end
end
if per_it_display==-1
    per_it_display = sprintf('%02d:%05.2f per iteration',t_min,t_sec);
else
    per_it_display = sprintf('%.1f iterations per second',1/t_per_it);
end

% Determine time elapsed str (skipping days/hours if 0).
if t_min_total<60
    elapsed_str = sprintf('%02d:%02d',t_min_total,t_sec_total);
else
    t_hour_total = floor(t_min_total/60);
    t_min_total = t_min_total-60*t_hour_total;
    if t_hour_total<24
        elapsed_str = sprintf('%02d:%02d:%02d',t_hour_total,t_min_total,t_sec_total);
    else
        t_day_total = floor(t_hour_total/24);
        t_hour_total = t_hour_total-24*t_day_total;
        elapsed_str = sprintf('%d/%02d:%02d:%02d',t_day_total,t_hour_total,t_min_total,t_sec_total);
    end
end

out_    = cell(1,3);
out_{1} = sprintf('%05.1f%% done after a total time of %s.', 100*n/N_total,elapsed_str);
out_{2} = sprintf('estimated time of completion: %s',        ETA_);
out_{3} = sprintf('(%s, %d iterations left)',                per_it_display,N_total-n);
if nargout==0
    clc
    fprintf('%s\n',out_{:});
end
if nargout>=1,out = out_;end
if nargout==2
    out_struct = struct(...
        'p_done',sprintf('%05.1f%%',100*n/N_total),...
        'elapsed_str',elapsed_str,...
        'ETA',ETA,...
        'ETA_date',ETA_date,...
        't_min',t_min,...
        't_sec',t_sec,...
        'it_per_sec',1/t_per_it);
end
end
function tf=ifversion(test,Rxxxxab,Oct_flag,Oct_test,Oct_ver)
%Determine if the current version satisfies a version restriction
%
% To keep the function fast, no input checking is done. This function returns a NaN if a release
% name is used that is not in the dictionary.
%
% Syntax:
%   tf = ifversion(test,Rxxxxab)
%   tf = ifversion(test,Rxxxxab,'Octave',test_for_Octave,v_Octave)
%
% Input/output arguments:
% tf:
%   If the current version satisfies the test this returns true. This works similar to verLessThan.
% Rxxxxab:
%   A char array containing a release description (e.g. 'R13', 'R14SP2' or 'R2019a') or the numeric
%   version (e.g. 6.5, 7, or 9.6). Note that 9.10 is interpreted as 9.1 when using numeric input.
% test:
%   A char array containing a logical test. The interpretation of this is equivalent to
%   eval([current test Rxxxxab]). For examples, see below.
%
% Examples:
% ifversion('>=','R2009a') returns true when run on R2009a or later
% ifversion('<','R2016a') returns true when run on R2015b or older
% ifversion('==','R2018a') returns true only when run on R2018a
% ifversion('==',23.02) returns true only when run on R2023b
% ifversion('<',0,'Octave','>',0) returns true only on Octave
% ifversion('<',0,'Octave','>=',6) returns true only on Octave 6 and higher
% ifversion('==',9.10) returns true only when run on R2016b (v9.1) not on R2021a (9.10).
%
% The conversion is based on a manual list and therefore needs to be updated manually, so it might
% not be complete. Although it should be possible to load the list from Wikipedia, this is not
% implemented.
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 1.2.1.1                                                       |%
%|  Date:    2023-10-20                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). You can see the full test matrix below.
% Compatibility considerations:
% - This is expected to work on all releases.

if nargin<2 || nargout>1,error('incorrect number of input/output arguments'),end

% The decimal of the version numbers are padded with a 0 to make sure v7.10 is larger than v7.9.
% This does mean that any numeric version input needs to be adapted. multiply by 100 and round to
% remove the potential for float rounding errors.
% Store in persistent for fast recall (don't use getpref, as that is slower than generating the
% variables and makes updating this function harder).
persistent  v_num v_dict octave
if isempty(v_num)
    % Test if Octave is used instead of Matlab.
    octave = exist('OCTAVE_VERSION', 'builtin');
    
    % Get current version number. This code was suggested by Jan on this thread:
    % https://mathworks.com/matlabcentral/answers/1671199#comment_2040389
    v_num = [100, 1] * sscanf(version, '%d.%d', 2);
    
    % Get dictionary to use for ismember.
    v_dict = {...
        'R13' 605;'R13SP1' 605;'R13SP2' 605;'R14' 700;'R14SP1' 700;'R14SP2' 700;
        'R14SP3' 701;'R2006a' 702;'R2006b' 703;'R2007a' 704;'R2007b' 705;
        'R2008a' 706;'R2008b' 707;'R2009a' 708;'R2009b' 709;'R2010a' 710;
        'R2010b' 711;'R2011a' 712;'R2011b' 713;'R2012a' 714;'R2012b' 800;
        'R2013a' 801;'R2013b' 802;'R2014a' 803;'R2014b' 804;'R2015a' 805;
        'R2015b' 806;'R2016a' 900;'R2016b' 901;'R2017a' 902;'R2017b' 903;
        'R2018a' 904;'R2018b' 905;'R2019a' 906;'R2019b' 907;'R2020a' 908;
        'R2020b' 909;'R2021a' 910;'R2021b' 911;'R2022a' 912;'R2022b' 913;
        'R2023a' 914;'R2023b' 2302};
end

if octave
    if nargin==2
        warning('HJW:ifversion:NoOctaveTest',...
            ['No version test for Octave was provided.',char(10),...
            'This function might return an unexpected outcome.']) %#ok<CHARTEN>
        if isnumeric(Rxxxxab)
            v = 0.1*Rxxxxab+0.9*fixeps(Rxxxxab);v = round(100*v);
        else
            L = ismember(v_dict(:,1),Rxxxxab);
            if sum(L)~=1
                warning('HJW:ifversion:NotInDict',...
                    'The requested version is not in the hard-coded list.')
                tf = NaN;return
            else
                v = v_dict{L,2};
            end
        end
    elseif nargin==4
        % Undocumented shorthand syntax: skip the 'Octave' argument.
        [test,v] = deal(Oct_flag,Oct_test);
        % Convert 4.1 to 401.
        v = 0.1*v+0.9*fixeps(v);v = round(100*v);
    else
        [test,v] = deal(Oct_test,Oct_ver);
        % Convert 4.1 to 401.
        v = 0.1*v+0.9*fixeps(v);v = round(100*v);
    end
else
    % Convert R notation to numeric and convert 9.1 to 901.
    if isnumeric(Rxxxxab)
        % Note that this can't distinguish between 9.1 and 9.10, and will the choose the former.
        v = fixeps(Rxxxxab*100);if mod(v,10)==0,v = fixeps(Rxxxxab)*100+mod(Rxxxxab,1)*10;end
    else
        L = ismember(v_dict(:,1),Rxxxxab);
        if sum(L)~=1
            warning('HJW:ifversion:NotInDict',...
                'The requested version is not in the hard-coded list.')
            tf = NaN;return
        else
            v = v_dict{L,2};
        end
    end
end
switch test
    case '==', tf = v_num == v;
    case '<' , tf = v_num <  v;
    case '<=', tf = v_num <= v;
    case '>' , tf = v_num >  v;
    case '>=', tf = v_num >= v;
end
end
function val=fixeps(val)
% Round slightly up to prevent rounding errors using fix().
val = fix(val+eps*1e3);
end

