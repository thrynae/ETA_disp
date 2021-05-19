function out=ETA_disp(N_total,n,then,per_it_display,scale_estimate)
%Display the estimated time of completion of a looping process.
%
% ETA_disp(N_total,n,then)
% ETA_disp(N_total,n,then,per_it_display)
% ETA_disp(N_total,n,then,per_it_display,scale_estimate)
% msg=ETA_disp(___)
%
%input:
%  N_total: total number of iterations
%  n: current count of iterations (will be adjusted up to eps if lower than that)
%  then: the output of the now function when the loop began
%  per_it_display: (optional, default is 0) 0 for automatic, -1 to display time per iteration on
%                             the last line, or 1 to print iterations per seconds on the last line.
%  scale_estimate: the time remaining is multiplied by this factor (this doesn't change the rate
%                                                         calculated over the processed iterations)
%
%output:
%  if called without output arguments, the message is printed to the command window, otherwise the
%  message is returned in a cell array.
%  The time remaining is capped at 10^10 seconds.
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
%       if mod(n,10^4)==0 %Avoid spending too much time printing the ETA.
%           ETA_disp(N,n,then,1)%Print iterations/s instead of s/iteration.
%       end
%   end
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 1.2.0                                                         |%
%|  Date:    2021-05-19                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - This is expected to work on all releases.

persistent legacy
if isempty(legacy)
    % The addtodate function was probably introduced in ML7.0, but the 'second' option was only
    % implemented later, at the latest in R2010a.
    legacy.addtodate=ifversion('<','R2010a','Octave','<',4);
end
if nargin<4 || numel(per_it_display)~=1,per_it_display=0;end
if nargin<5 || numel(scale_estimate)~=1,scale_estimate=1;end
n=max(eps,n); % Ensure n is non-zero, which would cause problems later.
t=(now-then)*(24*60*60); % This is equivalent to t=toc(h_tic);
t_min_total=floor(t/60);t_sec_total=round(t-60*t_min_total); % Calculate time elapsed.
f=n/N_total;t_r=(t-f*t)/f; % Calculation fraction done and time remaining.
if nargin>=5,t_r=t_r*scale_estimate;end % Scale the remainder estimate.
t_r=min(10^10,round(t_r)); % Round and cap to 10^10 seconds.
if legacy.addtodate
    ETA=datestr(now+(t_r/(24*60*60)),'HH:MM:SS');
else
    ETA=datestr(addtodate(now,t_r,'second'),'HH:MM:SS');
end
t_per_it=t/n; % This is in seconds per iteration.
if abs(per_it_display)~=1 % Revert to dynamic if the input is invalid.
    if t_per_it>0.5
        per_it_display=-1; % Display time per iteration.
    else
        per_it_display= 1; % Display iterations per second.
    end
end
if per_it_display==-1
    t_min=floor(t_per_it/60);t_sec=t_per_it-60*t_min;
    per_it_display=sprintf('%02d:%05.2f per iteration',t_min,t_sec);
else
    per_it_display=sprintf('%.1f iterations per second',1/t_per_it);
end

out=cell(1,3);
out{1}=sprintf('%05.1f%% done after a total time of %02d:%02d.',...
    100*n/N_total,t_min_total,t_sec_total);
out{2}=sprintf('estimated time of completion: %s',...
    ETA);
out{3}=sprintf('(%s, %d iterations left)',...
    per_it_display,N_total-n);
if nargout==0
    clc
    fprintf('%s\n',out{:});
    clear
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
% ifversion('==',9.9) returns true only when run on R2020b
% ifversion('<',0,'Octave','>',0) returns true only on Octave
% ifversion('<',0,'Octave','>=',6) returns true only on Octave 6 and higher
%
% The conversion is based on a manual list and therefore needs to be updated manually, so it might
% not be complete. Although it should be possible to load the list from Wikipedia, this is not
% implemented.
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 1.0.6                                                         |%
%|  Date:    2021-03-11                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - This is expected to work on all releases.

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
    ii=strfind(v_num,'.');if numel(ii)~=1,v_num(ii(2):end)='';ii=ii(1);end
    v_num=[str2double(v_num(1:(ii-1))) str2double(v_num((ii+1):end))];
    v_num=v_num(1)+v_num(2)/100;v_num=round(100*v_num);
    
    %get dictionary to use for ismember
    v_dict={...
        'R13' 605;'R13SP1' 605;'R13SP2' 605;'R14' 700;'R14SP1' 700;'R14SP2' 700;
        'R14SP3' 701;'R2006a' 702;'R2006b' 703;'R2007a' 704;'R2007b' 705;
        'R2008a' 706;'R2008b' 707;'R2009a' 708;'R2009b' 709;'R2010a' 710;
        'R2010b' 711;'R2011a' 712;'R2011b' 713;'R2012a' 714;'R2012b' 800;
        'R2013a' 801;'R2013b' 802;'R2014a' 803;'R2014b' 804;'R2015a' 805;
        'R2015b' 806;'R2016a' 900;'R2016b' 901;'R2017a' 902;'R2017b' 903;
        'R2018a' 904;'R2018b' 905;'R2019a' 906;'R2019b' 907;'R2020a' 908;
        'R2020b' 909;'R2021a' 910};
end

if octave
    if nargin==2
        warning('HJW:ifversion:NoOctaveTest',...
            ['No version test for Octave was provided.',char(10),...
            'This function might return an unexpected outcome.']) %#ok<CHARTEN>
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
    elseif nargin==4
        % Undocumented shorthand syntax: skip the 'Octave' argument.
        [test,v]=deal(Oct_flag,Oct_test);
        % Convert 4.1 to 401.
        v=0.1*v+0.9*fix(v);v=round(100*v);
    else
        [test,v]=deal(Oct_test,Oct_ver);
        % Convert 4.1 to 401.
        v=0.1*v+0.9*fix(v);v=round(100*v);
    end
else
    % Convert R notation to numeric and convert 9.1 to 901.
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
    case '==', tf= v_num == v;
    case '<' , tf= v_num <  v;
    case '<=', tf= v_num <= v;
    case '>' , tf= v_num >  v;
    case '>=', tf= v_num >= v;
end
end