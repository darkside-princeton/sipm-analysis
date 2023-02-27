# Quick Start

The analysis workflow is divided in two stages. The first stage is pre-processing of digitized waveforms, and the second stage is the higher-level analyses on the pre-processed data from the first stage. For more details, please read issue [#3](https://github.com/darkside-princeton/sipm-analysis/issues/3).

Data pre-processing is performed with the Python scripts under ``/sipm/exe/``. One could execute the scripts directly, or one could also submit jobs to the cluster to process a collection of data sets in parallel. The job submission is taken care of in ``jupyter/submit.ipynb``. Execute the cell for the script you want to run. Each of the scripts is used for the following tasks:
- ``analysis.py``: A minimal working example
- ``laser_pulse.py``: SiPM calibration
- ``scintillation_pulse.py``: LAr scintillation pulse analysis
- ``laser_waveform.py``: SiPM template analysis (under development)
- ``scintillation_waveform.py``: Triplet lifetime analysis (under development)

The script will produce an HDF5 file containing the processed data, such as amplitude, charge, and baseline rms. One can write one's own code to perform high-level analysis. There are also existing Jupyter Notebooks under ``jupyter/``:
- ``io.ipynb``: Demonstration of h5 I/O. Takes as input the files produced from ``analysis.py``.
- ``calibration.ipynb``: SiPM calibration. Takes as input the files produced from ``laser_pulse.py``.
- ``spectrum.ipynb``: Light yield measurement. Takes as input the files produced from ``scintillation_pulse.py``.
- ``triplet_lifetime.ipynb``: As the name suggests. Takes as input the files produced from ``laser_waveform.py`` and ``scintillation_waveform.py`` (under development).
