# physioData_preprocessing

getting started:
- setup python 2 evironment, for example using anaconda and install all necessary packages and dependencies as listed in the first lines of the scripts
- run tk_physionoise.py from your conda env to clean physio data
- run retroicor_progressbar_tk.py from your conda env to filter resting state time series for physio noise

input needed:
- physio data (txt): dummy_data/physio_data_resting_state.txt
- resting state time series (4D nii): ask Claudia/Esther for Sensemap2 data

tk_physionoise.py - settings:\
keep default settings

tk_physionoise.py - what the script does:
- adapted from Kelley et al., 2008, Plos One, https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0001751
- tk_physionoise.py = GUI (includes JS_PhysioNoise_pyLog.py)
- actual script = JS_PhysioNoise_pyLog.py
- preprocess physio data 
- extracts respiration and pulse phase values
- generates single respiration and pulse files containing phase values
- generated files include one column containing all phase values for the time series

tk_physionoise.py - output:
- RW = Downsampled respiratory spline waveform; 
- RVT = Respiratory volume over time based on peakdet peaks; 
- RRT = Respiratory rate over time; 
- CPttl = Cardiac peak based on the TTL pulse from the scanner; 
- CPd3 = Cardiac peak based on trough of the third derivative; 
- CPd3R = Cardiac peak R-wave estimate based on small peak of third derivative preceding the CPd3;
- CVT = Cardiac volume over time based on peakdet peaks; 
- CRTttl = Cardiac rate over time based on the TTL pulse from the scanner; 
- CRTd3 = Cardiac rate over time based on CPd3;
- CRTd3R = Cardiac rate over time based on CPd3R).

retroicor_progressbar_tk - settings:
- acq order: interleaved2 (equal # of slices → 2nd slice first / even slices first)
- TR: 2.0
- Freq: 200
- additional delay: 0.050s (pulse offset for measurement at the finger)
- IMPORTANT: Pulse was measured at the fingertip. The pulse wave from the heart arrives earlier in the brain than in the finger. Respiration on the other hand directly biases magnetic field and therefore appears immediately in the MR signal.
- find pulse offset by calculating residuals for different offsets: retroicor_tk_multidelay.py (generates 4D nii for different offsets)

retroicor_progressbar_tk.py - what the script does:
- includes js_validateentry.py
- cleans resting-state time series from physio noise
- use phase values (0 - 2Pi) to generate 2nd grade Fourier series: sin(x) + cos(x) + Sin(2x) + Cos(2x) for every single slice
- results go into GLM and add one column filled with “1” as average 
- residuals are taken as cleaned / filtered data
- additional pulse offset implemented, because pulse was measured at the fingertip with pulse oxi, pulse wave arrives at different times in the brain and the finger, shift the pulse phase results to minimize residuals in filtered epi data (until shifting is too much and variance in the data becomes greater again)

retroicor_progressbar_tk - generated graphics:
- look at Bodeplot and image 3 
- bodeplot shows crosstalk between pulse and breathing, e.g. heavy pulsation of abdominal aorta also visible in respiration signal obtained with breathing belt, particularly in very thin people; if pulse finger clip lies on the belly and cable is moved by mere abdominal breathing respiration shows up in pulse signal 
- in image 3 control whether maxima were set correctly









