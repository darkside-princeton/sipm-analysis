# How to run the code

The analysis workflow is divided into two stages. The first stage is pre-processing of digitized waveforms, and the second stage is the higher-level analyses on the pre-processed data from the first stage. For more details, please read issue [#3](https://github.com/darkside-princeton/sipm-analysis/issues/3).

Data pre-processing is performed with the Python scripts under ``/sipm/exe/``. One could execute the scripts directly or submit jobs to the cluster to process a collection of data sets in parallel. The job submission is taken care of in ``jupyter/submit.ipynb``, where one can execute the cell for the script one wants to run. Each of the scripts is used for the following task:
- ``analysis.py``: A minimal working example
- ``laser_pulse.py``: SiPM calibration
- ``scintillation_pulse.py``: LAr scintillation pulse analysis
- ``laser_waveform.py``: SiPM template analysis (under development)
- ``scintillation_waveform.py``: Triplet lifetime analysis (under development)

Each script will produce an HDF5 file under ``/scratch/gpfs/your-netid/results/`` containing the processed data, including amplitude, charge, baseline rms, etc. One could write one's own high-level analysis code or use the existing Jupyter Notebooks under ``jupyter/``:
- ``io.ipynb``: Demonstration of h5 I/O. Takes as input the files produced from ``analysis.py``.
- ``calibration.ipynb``: SiPM calibration. Takes as input the files produced from ``laser_pulse.py``.
- ``spectrum.ipynb``: Light yield measurement. Takes as input the files produced from ``scintillation_pulse.py``.
- ``triplet_lifetime.ipynb``: As the name suggests. Takes as input the files produced from ``laser_waveform.py`` and ``scintillation_waveform.py`` (under development).
