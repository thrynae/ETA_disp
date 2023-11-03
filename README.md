# ifversion documentation
[![View ifversion on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/69138-ifversion)
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=thrynae/ifversion)

Table of contents

- Description section:
- - [Description](#description)
- Matlab/Octave section:
- - [Syntax](#syntax)
- - [Output arguments](#output-arguments)
- - [Input arguments](#input-arguments)
- - [Examples](#examples)
- - [Compatibility, version info, and licence](#compatibility-version-info-and-licence)

## Description
This function will print a short message to the command window, giving you an idea when a process might finish. It displays the current percentage done, the runtime until now, the estimated time of completion, the time per iteration (or iterations per second), and the number of iterations left.

Alternatively you can capture the text before it is printed and use the text in your own application, e.g. in a GUI.

## Matlab/Octave

### Syntax

    ETA_disp(N_total,n,then)
    ETA_disp(N_total,n,then,per_it_display)
    ETA_disp(N_total,n,then,per_it_display,scale_estimate)
    msg = ETA_disp(___)
    [msg,strct] = ETA_disp(___)

### Output arguments

|Argument|Description|
|---|---|
|msg||
|strct||

### Input arguments

|Argument|Description|
|---|---|
|N_total|The total number of iterations.|
|n|The current count of iterations (this will be adjusted up to `eps` if 0).|
|then|This should be the value the `now` function returned before the loop began.|
|per_it_display| (optional, default is 0) if set to -1, the function displays the time per iteration on the last line. Setting to 1 will print iterations per seconds on the last line. With 0 the choice will be made automatically.|
|scale_estimate|The time remaining is multiplied by this factor (this doesn't change the rate calculated over the processed iterations). This allows skewing the estimate in cases where the earlier iterations are expected to be faster.|

### Examples

    N = 10;then = now; %#ok<TNOW1>
    for n=1:N
        pause((1+rand)*2/N)
        ETA_disp(N,n,then)
    end
 
    N = 10^6;then = now; %#ok<TNOW1>
    for n=1:N
        pause((1+rand)*2/N)
        if mod(n,10^4)==0 % Avoid spending too much time printing the ETA.
            ETA_disp(N,n,then,1) % Print iterations/sec instead of sec/iteration.
        end
    end


### Compatibility, version info, and licence
Compatibility considerations:
- This is expected to work on all releases.

|Test suite result|Windows|Linux|MacOS|
|---|---|---|---|
|Matlab R2023b|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it>Monterey : Pass</it>|
|Matlab R2023a|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Matlab R2022b|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it>Monterey : Pass</it>|
|Matlab R2022a|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Matlab R2021b|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it>Monterey : Pass</it>|
|Matlab R2021a|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Matlab R2020b|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it>Monterey : Pass</it>|
|Matlab R2020a|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Matlab R2019b|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it>Monterey : Pass</it>|
|Matlab R2019a|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Matlab R2018a|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it></it>|
|Matlab R2017b|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it>Monterey : Pass</it>|
|Matlab R2016b|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it>Monterey : Pass</it>|
|Matlab R2015a|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it></it>|
|Matlab R2013b|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Matlab R2007b|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Matlab 6.5 (R13)|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Octave 8.2.0|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Octave 7.2.0|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Octave 6.2.0|<it>W11 : Pass</it>|<it>raspbian_11 : Pass</it>|<it>Catalina : Pass</it>|
|Octave 5.2.0|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Octave 4.4.1|<it>W11 : Pass</it>|<it></it>|<it>Catalina : Pass</it>|

    Version: 2.0.0
    Date:    2023-11-03
    Author:  H.J. Wisselink
    Licence: CC by-nc-sa 4.0 ( https://creativecommons.org/licenses/by-nc-sa/4.0 )
    Email = 'h_j_wisselink*alumnus_utwente_nl';
    Real_email = regexprep(Email,{'*','_'},{'@','.'})

### Test suite

The tester is included so you can test if your own modifications would introduce any bugs. These tests form the basis for the compatibility table above. Note that functions may be different between the tester version and the normal function. Make sure to apply any modifications to both.
