# PipeSound Library     

&nbsp; 

The PipeSound Library is a collection of Matlab scripts that includes some of the software blocks of the audio synthesiser used to generate the [PipeSound](https://doi.org/10.5281/zenodo.7615371) synthetic dataset. The repo contains the following folders:

- An ***acoustic model*** to calculate the dispersive behaviour of a waveguide composed of a cylindrical elastic shield filled with inviscid liquid,
- A ***noise filter*** to improve the signal-to-noise ratio of real recordings,
- A ***set of examples*** for the acoustic model and the noise filter. A snippet of the [PipeSound](https://doi.org/10.5281/zenodo.7615371) soundbank is also included to provide a small collection of real recordings to feed the noise filter.

The library includes the following items:

- ``nameSpace.m``: .m file containing the definitions of names and paths required by the library,  
- ``spaceSettings.m``: .m file containing the default values,  
 ``pathManager.m``: .m file adding the necessary paths for the Matlab environment,  
- ``AcousticModel``: a folder containing the scripts required for the acoustic model,  
- ``NoiseFilter``: a folder containing the scripts required for the noise filter,  
- ``Examples``: a folder containing the examples and the soundbank snippet,  
- ``Utilities``: a folder containing some general-purpose utility functions.  

&nbsp;   
>**Note**
>To use the library, download the repository and set ``ns.repository.path`` in ``nameSpace.m`` to the path of the ``Examples`` folder. No path needs to be added to Matlab before running the scripts.


&nbsp;    
&nbsp;  
## Main files   

A brief description of the main scripts is reported below. A comprehensive description is included at the top of each file. Executable scripts are marked with *.   

&nbsp;   
``AcousticModel``:

- ``kzSolver.m*``: a numerical solver for the calculation of the axial wavenumbers *k<sub>z</sub>*,
- ``modeExtractor.m*``: a mode extractor for the results obtained from ``kzSolver.m``,
- ``characteristicMatrix.m``: characteristic matrix of the waveguide used by ``kzSolver.m`` to calculate the axial wavenumbers *k<sub>z</sub>*.

&nbsp;  
``NoiseFilter``:

- ``getNoiseFrequencyMask.m``: calculates the noise mask for the noise filter,
- ``noiseFilter.m``: a spectral filter to improve the signal-to-noise ratio of real recordings,
- ``filterExample.m*``: a simple script using ``noiseFilter.m`` to filter real recordings from the soundbank snippet,
- ``signalNoiseRatioTestUtility.m*``: a script to measure the filter performance in terms of noise and signal attenuation and signal-to-noise ratio gain.

&nbsp;  
``Examples``:

Examples are provided for the acoustic model and the noise filter. The files for the required settings can be found in the related folders.


&nbsp;    
&nbsp; 

## Contacts   

Further requests: <pipesoundlibrary@gmail.com>





