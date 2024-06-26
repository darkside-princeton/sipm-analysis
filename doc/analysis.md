# How to run the code

The analysis workflow is divided into two stages. The first stage is pre-processing of digitized waveforms, and the second stage is the higher-level analyses on the pre-processed data from the first stage. For more details, please read issue [#3](https://github.com/darkside-princeton/sipm-analysis/issues/3).

# Pre-processing

Data pre-processing is performed with the Python scripts under ``sipm/exe/``. One could execute the scripts directly or submit jobs to the cluster to process a collection of data sets in parallel. Job submissions can be done using ``sipm/submit/submit.py``. This program takes as input a configuration file under ``sipm/config/``. For example,

```
python sipm/submit/submit.py -c sipm/config/analysis.yaml
```

| Config. File | Description | Higher-level analysis | Pre-processing |
| ------------ | ----------- | --------------------- | -------------- |
| sipm/config/analysis.yaml | A minimal example | jupyter/io.ipynb | sipm/exe/analysis.py |
| sipm/config/laser_pulse_*.yaml | Laser calibration | jupyter/calibration.ipynb | sipm/exe/laser_pulse.py |
| sipm/config/laser_waveform_*.yaml | SiPM single PE pulse shape | jupyter/triplet_lifetime.ipynb jupyter/deconvolution.ipynb | sipm/exe/laser_waveform.py |
| sipm/config/scintillation_pulse_*.yaml | Scintillation light from Cs-137 gammas | jupyter/spectrum.ipynb | sipm/exe/scintillation_pulse.py |
| sipm/config/scintillation_waveform_*.yaml | LAr scintillation pulse shape | jupyter/triplet_lifetime.ipynb | sipm/exe/scintillation_waveform.py |
| sipm/config/fft_*.yaml | FFTs on noise and single PE waveforms | jupyter/deconvolution.ipynb | sipm/exe/fft.py

Each script will produce HDF5 files under ``/scratch/gpfs/your-netid/results/`` containing the processed data, including amplitude, charge, baseline rms, etc.
The file names end with the name of the pre-processing script + .h5.

# Higher-level analysis

One could write one's own high-level analysis code or use the existing Jupyter Notebooks under ``jupyter/``:
- ``io.ipynb``: Demonstration of h5 I/O. Takes as input the files produced from ``analysis.py``.
- ``calibration.ipynb``: SiPM calibration. Takes as input the files produced from ``laser_pulse.py``.
- ``spectrum.ipynb``: Light yield measurement. Takes as input the files produced from ``scintillation_pulse.py``.
- ``triplet_lifetime.ipynb``: As the name suggests. Takes as input the files produced from ``laser_waveform.py`` and ``scintillation_waveform.py``.
- ``deconvolution.ipynb``: Testing different deconvolution methods and developing pulse finding algorithms. Work in progress.

In the future, we will implement Python classes for these analysis code.
