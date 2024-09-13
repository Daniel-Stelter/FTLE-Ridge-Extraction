This repository provides implementations used in the following submission:

Stelter, D., Wilde, T., Rössl, C., & Theisel, H. (2024). A Particle-Based Approach to Extract Dynamic 3D FTLE Ridge Geometry. Computer Graphics Forum.

These implementations contain the following methods:

- our proposed particle system,
- the particle system of Kindlmann et al.[^1],
- and the Marching Ridges method of Furst and Pizer[^2],

which are evaluated and compared in our work.

# Installation

For building the project, you must have installed the following prerequisites to your system:

- C++ compiler implementing the C++17 standard
- CMake (at least version 3.9)
- OpenMP library (if not included with the compiler)
- Eigen3 library

For downloading the code, please be aware that this project uses a submodule.
Thus, you can either clone this repository as follows ...
```
git clone --recurse-submodules https://github.com/Daniel-Stelter/FTLE-Ridge-Extraction.git
```
... or if you have the code already downloaded, use the following code inside the main project directory:
```
git submodule update --init --recursive
```
After this, you can follow the typical installation steps using the `CMakeLists.txt` file.
Assuming you are currently in the project directory:
```
mkdir build
cd build
cmake ..
make
```
This generates three executables which are explained in the next section.

**Tested systems:**
| OS | Compilers |
|----|-----------|
| Ubuntu 22.04 | Clang 14.0.0 and GCC 11.4.0 |
| Arch Linux 2024.08.01 | Clang 18.1.8 and GCC 14.2.1 |
| macOS Sonoma 14.6.1 | Clang 15.0.0 (Apple) |


# Execution

This project mainly proposes executions for extracting FTLE ridges from flow datasets.
Our project provides two build targets:

- `execute`
- `test-ridges`

The `execute` target is the main project which applies the extraction methods to flow datasets.
The `test-ridges` target performs extractions to a test example.
Both are explained in the following.


## Ridge extraction from flow datasets

The `execute` target applies one of the three methods to a dataset.
It takes one argument which is a path to a setup file which defines all parameters for the experiment.
We already deliver a multitude of setups which were used for the submitted work.
These can be found in the `setups` directory in our project.
One can either apply the extraction methods to implemented analytical datasets, or to flows saved as an Amira file (more details below).


### Quick start

After successfully compiling the code, one can use a readily available setup file.
For example, given one is located in a build directory inside the main project, the execution could look like:

```
./execute ../setups/Stelter/ABC3D.json
```

This will generate all results inside a directory which is created inside the current working directory.

In order to use Amira flows, one has to download the respective files.
For the existing setup files, you can find the respective data here:

- [Half Cylinder](https://cgl.ethz.ch/Downloads/Data/ScientificData/halfcylinder3d-Re320_am.zip)
- [Cloud-topped Boundary Layer](https://cgl.ethz.ch/Downloads/Data/ScientificData/ctbl3d_am.zip)

**Note:** Please make sure that the location of your extracted files and the paths specified in the setup file match.
The paths in the file are relative w.r.t. the current working directory.


### Create / change setup files

We use JSON files for the project setup.
There is a multitude of parameters which have to be specified.
If not specifically declared otherwise, it is not permissible to leave out any fields.

| Field | Description | Format/Values |
|-------|-------------|---------------|
| `dataset` | Name of a dataset. Can either be an implemented dataset or an Amira dataset. These names **must** end with either `2D` or `3D`, and in case of Amira datasets they **must** begin with `Amira-`. | `ABC3D`, `DoubleGyre2D`, `DoubleGyre3D`, `ForcedDuffing2D`, `Tornado3D`, `Amira-...2D`, `Amira-...3D` |
| `amira_files` | Only required for Amira datasets. This list has to contain all save files **in the correct order**. | String[] |
| `method` | Method for extracting the ridges. Please note: `Marching Ridges` is only applicable to 3D datasets. | `Stelter`, `Kindlmann`, `Marching Ridges` |
| `save_dir` | Save directory for the output files. Supports '~' as symbol for the home directory. | String |
| `domain` | Domain used for the ridge extraction. `null` resolves to the defined (or for analytical datasets, internally selected) domain. | `[[xmin,xmax],[ymin,ymax]{,[zmin,zmax]}]`, `null` |
| `t` | Start time for the experiment. | Float |
| `tau` | Maximum integration time for the experiment. | Float |
| `steps` | Integration time steps for the experiment. | Integer |
| `space_res` | Resolution for the FTLE field. | `[x,y{,z}]` |
| `min_ridge_strength` | Threshold for ridge strength (i.e., sharpness of ridges) for the extraction. | Float |
Further parameters for 2D datasets:
| `img_res` | Resolution for the output images. | `[u,v]` |
| `rendered_domain` | Domain which should be rendered. `null` resolves to the value of `domain`. | `[[xmin,xmax],[ymin,ymax]]`, `null` |
Further paramters for Stelter:
| `voxel_res` | Resolution for the voxel grid (for internally saving particles). | `[x,y{,z}]` |
| `par_init_res` | Resolution for particle initialization. | `[x,y{,z}]` |
| `full_init_period` | Number of steps after which new particles are initialized. | Integer |
| `elliptic_dis_fac` | Number of steps after which new particles are initialized. | Integer |
Further paramters for Kindlmann:
| `voxel_res` | Resolution for the voxel grid (for internally saving particles). | `[x,y{,z}]` |
| `par_init_res` | Resolution for particle initialization. | `[x,y{,z}]` |
| `sigma` | Radius for the neighborhood space of a particle. | Float |
| `numeric_diff_delta` | Spatial distance for numerical computation of gradients and Hessians. | Float |


### Notes on Amira datasets

One can execute FTLE ridge extractions to arbitrary Amira flow files / vector fields.
This includes time-dependent and time-independent vector fields.
However, there are some rules to follow:

- The `dataset` field in the setup file **must** have the prefix `Amira-`, and either the postfix `2D` or `3D`. Example: `Amira-MyVectorField2D`.
- One has to provide the `amira_files` field. Even if there is only one file, it has to be inside a list of Strings (e.g. `["/data/vector_field.am"]`).
- If the dataset consists of multiple files, they **must be provided in the correct order** in this list.


## Test executions

The traget `test-ridges` does not take any arguments.
It handles an analytical dataset with ridges which move in a circular shape.
This dataset is available for both 2D and 3D.
The execution as-is generares the following results:

- `ParSys-Circles2D-paths`: comparison of the constraining steps of our own and Kindlmann's particle systems
- `Stelter-Circles2D`: extraction of our method for the 2D version
- `Kindlmann-Circles2D`: extraction of the Kindlmann method for the 2D version
- `Stelter-Circles3D`: extraction of our method for the 3D version
- `Kindlmann-Circles3D`: extraction of the Kindlmann method for the 3D version
- `MarRidges-Circles3D`: extraction of the Marching Ridges method for the 3D version


# Obtaining the results

All results are generated in the ouput directory specified in the setup file.
For the setups delivered in our project this is `./results/<method>/<setup-name>`.
However, this can also be customized.
Inside this directory, multiple subdirectories are generated.
The `data-saves` directory contains data generated during the execution of the respective methods.
The `out` directory contains the output data of interest.
In case of 2D datasets, it contains image files of the basic FTLE field, and additionally images of the extracted ridges.
For 3D datasets, instead, it creates a forlder with objects which can be loaded in Blender.
To actually use these data, you need to use the following Blender scene:

[https://cloud.ovgu.de/s/5bMSmSi8ALsjJLm](https://cloud.ovgu.de/s/LkxNXrR9aQdWceg)

This scene contains the camera pose and light sources of our ABC flow example.
Using Blender 4.2.1, you can load data by following these steps:

1. Navigate to the `Scripting` tab and select the `load_data.py` script.
2. Change the `input_dir` path to your respective save file directory.
3. Depending on whether you load data of a particle system or Marching Ridges, comment / uncomment the respective lines:
   - `loadParSys` and `loadMarRidges` load an object of a single time step
   - `renderParSysTimeSeries` and `renderMarRidgesTimeSeries` automatically load, render and unload all objects inside the input directory
4. Optionally, you can change the size of the particles: Navigate to the `Geometry Nodes` tab and change the `Particle size` value. Please note: This also has an effect on all renderings of the time series rendering methods!


# References

[^1]: Kindlmann, G. L., Estépar, R. S. J., Smith, S. M., & Westin, C. F. (2009). Sampling and visualizing creases with scale-space particles. IEEE transactions on visualization and computer graphics, 15(6), 1415-1424.

[^2]: Furst, J. D., & Pizer, S. M. (2001). Marching ridges. Signal and Image Processing.
