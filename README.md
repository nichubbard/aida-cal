# aida-cal
Calibrate (gain match) AIDA from in-beam data

Based on M. Reese, et al., "Automatic instrinsic calibration of double-sided silicon detectors"
(https://doi.org/10.1016/j.nima.2015.01.032)

A python script to perfomr the intrinsic gain calibration of AIDA (DESPEC) via in-beam data.
Requires text data of light ions produced from the DESPEC Go4 code (or similar) of the format:

```
DSSD XStrip XAmplitude YStrip YAmplitude
```

Amplitude should be offset corrected (from a pulser walkthrough) and sign corrected 

There is no limit to the number of lines in the file - more is better

Suggested conditions are E > 1 MeV and single strip fired on both X and Y only
This should be selective for the light-ions.

## Usage
./aida-cal.py [options]

./aida-cal.py --help prints all the options.

A suggested workflow would be

```
./aida-cal.py --dssds=N
``` 

To generate the relative gains
Then analyse alpha background data with these gains (copy aida_gains.txt to
Go4/Configuration Files/AIDA/AIDA_gains.txt) and find the high energy 214Po Peak
Find the given energy of this peak and rerun

```
./aida-cal.py --use-gains --dssds=N --absolute=AAAA=7687
```

Where AAAA is the energy (keV) of the 214Po in the relatively matched file

This final aida_gains.txt will then contain absolute calibrated gains
for use in analysis


