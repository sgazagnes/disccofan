Usage
=====

.. role:: bold-smallcaps

:bold-smallcaps:`disccofan` can be utilized in various ways: through a configuration file, command-line options, with multi-threading or multiple MPI processes... Below are some basic examples of how to use :bold-smallcaps:`disccofan`:

1. **Example 1: Using the default `config.ini` file with a single process**:

   .. code-block:: bash

      ./disccofan 

2. **Example 2: Using a specific configuration file with a single process**:

   .. code-block:: bash

      ./disccofan -config myconfig.ini

3. **Example 3: Using the default `config.ini` file while overriding some parameters from the terminal**:

   .. code-block:: bash

      ./disccofan -input_image ilikebirds.png 

4. **Example 4: Using 4 MPI processes, 2 threads, and the default `config.ini` file**:

   .. code-block:: bash

      mpirun -n 4 ./disccofan -g 2,2,1 -threads 2

:bold-smallcaps:`disccofan` provides a wide range of parameters and options for customizing the processing of data sets. Below is a detailed description of the configuration parameters and command-line options.

Configuration File Parameters
-----------------------------

- **input_name**: (str) The name of the input file.
- **output_name**: (str) The name of the output file.
- **input_type**: (str) File type for the input.
- **output_type**: (str) File type for the output.
- **hdf5_dataset**: (str) The dataset name if using HDF5 files.
- **save_output**: (int) Whether to save output images (1 for yes, 0 for no).
- **tile_overlap**: (int) Indicates if tiles include overlap (1 for yes, 0 for no).
- **mpi_grid**: (list of int) Number of tiles horizontally, vertically, and in depth.
- **threads**: (int) Number of threads to use.
- **preprocessing**: (str) Preprocessing of the data.
- **image_options**: (str) Additional preprocessing options.
- **pixel_dim**: (list of float) Physical dimensions of pixels.
- **operation**: (tuple) Processing operations, including filters and parameters.
- **tree_type**: (str) Type of tree operation.
- **pixel_connectivity**: (int) Pixel connectivity (4 for 2D, 6 for 3D, 8 or 26 for including diagonal neighbors).
- **verbosity**: (str) Logging level.

Command-Line Options
--------------------

- **-config**: (str) Path to the configuration file.
- **-input_image**: (str) Specify the input image file.
- **-output_image**: (str) Specify the output image file.
- **-input_type**: (str) Override input file type.
- **-output_type**: (str) Override output file type.
- **-hdf5_dataset**: (str) Override HDF5 dataset name.
- **-image_options**: (str) Override image options.
- **-verbosity**, **-v**: (str) Override verbosity level.
- **-save_output**: (int) Override save output setting.
- **-tile_overlap**: (int) Override tile overlap setting.
- **-threads**, **-t**: (int) Override number of threads.
- **-preprocessing**: (str) Override preprocessing setting.
- **-grid**, **-g**: (list of int) Override MPI grid settings.
- **-operation**, **-o**: (tuple) Override processing operations.

Detailed Parameter Descriptions
-------------------------------

**input_name**: Specifies the dataset to process. For command-line overrides, use `-input_image` (e.g., `-input_image toto.png`). Command-line options override those in the config file.

**input_type**: Leave as 'auto' for common formats (e.g., PNG, PGM, TIF, JP(E)G, FITS, HDF5).

For HDF5 files, specify the dataset name using **hdf5_dataset**.

For large datasets split into multiple parts, ensure enough MPI processes are allocated. Use the format `toto-T0T.[format]` for `input_name`, where `[format]` is the image format. `-T0T` means the data set is split and each part should be read by independant processes. See one of the example in the documentation for further details.

**tile_overlap**: Only for pre-cut data sets. Set to 1 if tiles overlap by one pixel; otherwise, set to 0. For single-file datasets, **tile_overlap** should be left to 1.

**pixel_dim**: Represents physical pixel dimensions if morphological accuracy is required.

**pixel_connectivity**: Specifies neighboring pixel connectivity. Default is 4 in 2D and 6 in 3D. Set to 8 or 26 to include diagonal neighbors.

Output File Parameters
----------------------

The **output_name** parameter in the configuration file specifies the name of the output file. The final name depends on the operation and parameters used. You can also specify the output file using the **-output_image** command-line option (e.g., `-output_image toto_out.png`).

For HDF5 output files, the dataset name defaults to `disccofan`.

Pre-processing Parameters
-------------------------

You can set the **preprocessing** option in the configuration file or via command line to preprocess data. Available options include:

- **ubyte**: 8-bit scaling
- **ushort**: 16-bit scaling
- **uint**: 32-bit scaling
- **log**: Natural logarithm
- **log10**: Logarithm base 10
- **exp**: Exponential
- **sqrt**: Square root

Although :bold-smallcaps:`disccofan` primarily handles grayscale datasets, it also supports RGB (converted to luminance), LOFAR, or specific sequences. Set **image_options** to "grayscale" (default), "rgb", "lofar", or "sequence" as needed.

Operation Parameters
--------------------

The **operation** parameter is crucial for specifying which analysis you want to perform on your data. There are eight possible operations: `filter`, `extract`, `dmp`, `csl`, `pattern`, `pattern2D`, `tree`, and `check`. Each operation requires an attribute function to specify the morphological properties of interest. Below, we describe these attributes and the corresponding operations.

Attributes
~~~~~~~~~~

There are 18 attribute functions categorized into 5 groups. Each function is indexed for reference and can be used to analyze the morphological properties of objects in the image.

1. **Group #1: Connected Components Size**

   - **Size**: Determines the size of objects in the image. Index: 0.

2. **Group #2: Connected Components Dimensions**

   Provides detailed measurements of object dimensions:

   - **Area of Minimal Enclosing Rectangle**: Area of the smallest rectangle that completely encloses the object. Index: 1.
   - **Square of Diagonal of Minimal Enclosing Rectangle**: Square of the diagonal length of the minimal enclosing rectangle. Index: 2.
   - **X-Extent**: Maximum extent of the object along the x-axis. Index: 3.
   - **Y-Extent**: Maximum extent of the object along the y-axis. Index: 4.
   - **Z-Extent**: Maximum extent of the object along the z-axis. Index: 5.

3. **Group #3: Connected Components Center of Mass**

   Provides coordinates of the center of mass:

   - **X-Mean**: X-coordinate of the center of mass. Index: 6.
   - **Y-Mean**: Y-coordinate of the center of mass. Index: 7.
   - **Z-Mean**: Z-coordinate of the center of mass. Index: 8.

4. **Group #4: Connected Components Weighted Center of Mass**

   Provides weighted center of mass information:

   - **X-Wmean**: X-coordinate of the weighted center of mass. Index: 9.
   - **Y-Wmean**: Y-coordinate of the weighted center of mass. Index: 10.
   - **Z-Wmean**: Z-coordinate of the weighted center of mass. Index: 11.
   - **Total Flux**: Total flux of the object. Index: 12.

5. **Group #5: Connected Components Inertia Properties**

   Provides details about the object's inertia properties:

   - **Trace of the Inertia Matrix**: Trace of the inertia matrix. Index: 13.
   - **Elongation**: Ratio of the longest axis to the second longest axis. Index: 14.
   - **Flatness**: Ratio of the second longest axis to the third longest axis (only in 3D). Index: 15.
   - **Sparseness**: Sparseness or porosity of the object. Index: 16.
   - **Non-compactness**: Non-compactness or asymmetry of the object. Index: 17.


Operations details
~~~~~~~~~~~~~~~~~~


1. **Filtering and Extracting Features**: 

   - **filter** or **extract** apply morphological operations based on chosen attributes. Example syntax:

      .. code-block:: bash
         
         operation = (
            ("filter", 0, 100)
         )
          
   This filters (or extracts) objects based on their size. `0` represents the attribute function index, and `100` is the threshold.

2. **Differential Morphological Profile (DMP)**: 

   - **dmp** generates a profile based on a sequence of extract operations with varying thresholds. Example syntax:

      .. code-block:: bash
   
         operation = (
            ("dmp", 0, "lvec.txt")
         )

   Here, `0` sets the attribute function, and `lvec.txt` contains the thresholds for filtering.

3. **CSL Segmentation**: 

   - **csl** segments the image based on scale, contrast, and luminance. Example syntax:

      .. code-block:: bash
   
         operation = (
            ("csl", 0, "lvec.txt")
         )

   `0` sets the attribute function, and `lvec.txt` contains thresholds for CSL properties. Note: This operation is partially implemented.

4. **Pattern Spectrum**: 

   - **pattern** calculates the 1D pattern spectrum. Example syntax:

      .. code-block:: bash
   
         operation = (
            ("pattern", 0, "lvec.txt", 1)
         )

   `0` sets the attribute function, `lvec.txt` contains thresholds for defining scales, and `1` is the scale factor applied to values in `lvec.txt`.

5. **2D Pattern Spectrum**: 

   - **pattern2D** calculates the 2D pattern spectrum. Example syntax:

      .. code-block:: bash
   
         operation = (
            ("pattern2D", 17, "lvec2.txt", 1, 0, "lvec2.txt", 1)
         )

   Here lvec.txt contains thresholds for attribute function 0, lvec2.txt contains thresholds for attribute function 17, and the two 1s are scale factors applied to values in each file.

6. **Tree**: 

   - **tree** generates a text file with the tree structure and attributes. Example syntax:

      .. code-block:: bash
   
         operation = (
            ("tree", 0)
         )

   Note: Using the tree operation with MPI may produce different results compared to sequential processing. Ensure you understand MPI parallelization before using it for tree analysis.

7. **Check**: 

   - **check** replaces each pixel value with the corresponding attribute value. Example syntax:

      .. code-block:: bash
   
         operation = (
            ("check", 0)
         )

   This fills each pixel with the size (attribute function 0) of the corresponding object it belongs to.


Combining Several Operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can combine multiple operations in a single run to avoid the overhead of rebuilding the tree structure and attributes for each operation. To do this, set the **operation** parameter with a list of operations. For example:

   .. code-block:: bash
      
      operation = (
         ("filter", 0, 100),
         ("filter", 13, 5),
         ("extract", 5, 1000),
         ("dmp", 1, "tests/lvec.txt"),
         ("pattern", 1, "tests/lvec.txt"),
         ("tree", 1)
      )

Operations will be executed sequentially. If operations use attribute functions from the same group, the computation will be faster since the same tree is reused.

Parallelization Options
-----------------------

You can enhance the performance of :bold-smallcaps:`disccofan` using threads or MPI. Threads are suitable for datasets with a narrow dynamic range (up to 16 bits). Set the number of threads using **threads** (or **-threads**).

For MPI parallelization, use `mpirun` or `mpiexec`. Define data distribution with the **mpi_grid** option (or **-grid**). For example, setting **mpi_grid** to `[2,2,1]` indicates 4 MPI processes. Ensure that the total number of MPI processes matches the grid size.

Tree Structure
--------------

By default, :bold-smallcaps:`disccofan` constructs a max-tree, emphasizing bright objects in the image. To focus on dark objects (with lower gray values), set the **tree_type** option to `min` to build a min-tree.

Verbosity
---------

Adjust the logging level using the verbosity option. Available levels include `off`, `debug`, `info`, `warn`, and `error`.

For additional help, refer to the examples section to see practical applications of these options.



