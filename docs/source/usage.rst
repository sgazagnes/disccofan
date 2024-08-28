Usage
=====

.. role:: bold-smallcaps


:bold-smallcaps:`disccofan` can work in many different ways, using a config file, command-line options, with multi-threading, or several MPI processes.
One easy way to use :bold-smallcaps:`disccofan` is through the config file. 

Here are some basic example usages of disccofan:

1. **Example 1: Using the default config.ini file and a single process**:

   .. code-block:: bash

      ./disccofan 

2. **Example 2: Using a specific config file and a single process**:

   .. code-block:: bash

      ./disccofan -config myconfig.ini


3. **Example 3: Using the default config.ini file and overriding some parameters from the terminal**:

   .. code-block:: bash

      ./disccofan -input_image ilikebirds.png 

4. **Example 4: Using 4 mpi processes, 2 threads, and the default config.ini file**:

   .. code-block:: bash

      mpirun -n 4 ./disccofan -g 2,2,1 -threads 2


:bold-smallcaps:`disccofan` offers an extensive range of parameters and options for fine-tuning the processing of the data sets. Below is a detailed description of the configuration parameters and command-line options.

Configuration File Parameters
-----------------------------		

- **input_name**: (str) The name of the input file.
- **output_name**: (str) The name of the output file.
- **input_type**: (str) File type for the input. 
- **output_type**: (str) File type for the output. 
- **hdf5_dataset**: (str) The dataset name if using HDF5 files.
- **save_output**: (int) Whether to save output images. 
- **tile_overlap**: (int) Whether tiles include overlap. 
- **mpi_grid**: (list of int) Number of tiles horizontally, vertically, and in depth.
- **threads**: (int) Number of threads to use.
- **preprocessing**: (str) preprocessing of the data.
- **image_options**: (str) Additional preprocessing options. 
- **pixel_dim**: (list of float) Physical dimensions of pixels.
- **operation**: (tuple) Processing operations, including filters and parameters.
- **tree_type**: (str) Type of tree operation. 
- **pixel_connectivity**: (int) Pixel connectivity. 
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
- **-verbosity, -v**: (str) Override verbosity level.
- **-save_output**: (int) Override save output setting.
- **-tile_overlap**: (int) Override tile overlap setting.
- **-threads, -t**: (int) Override number of threads.
- **-preprocessing**: (str) Override preprocessing setting.
- **-grid, -g**: (list of int) Override MPI grid settings.
- **-operation, -o**: (tuple) Override processing operations.


This is a lot to take, so we go into detailing each parameter hereafter.

Input File Parameters
---------------------

In the configuration file, **input_name** specifies the name of the dataset to process. You can also specify the dataset using the **-input_image** command-line option (e.g., `-input_image toto.png`). Note that command-line options always override those set in the config file.

:bold-smallcaps:`disccofan` automatically determines the data type and properties if the data format is common (e.g., PNG, PGM, TIF, JP(E)G, FITS, or HDF5). The **input_type** option should be left as 'auto'.

For HDF5 files, specify the dataset name with the **hdf5_dataset** parameter.

When dealing with large datasets divided into multiple parts to form a complete 2D or 3D image, ensure you have enough MPI processes for each part. The **input_name** should follow the format `toto-T0T.[format]`, where `[format]` is the image format. The `-T0T` suffix indicates pre-cut sections of the dataset.

The **tile_overlap** parameter is crucial in such cases. Set **tile_overlap** to 1 if tiles overlap by one pixel; otherwise, set it to 0. This ensures correct tile processing. For single-file datasets, you can leave **tile_overlap** as is.

The **pixel_dim** option reflects the physical size of pixels if needed for morphological accuracy.

**pixel_connectivity** specifies pixel neighboring connectivity. Default is 4 in 2D and 6 in 3D. Set to 8 or 26 to include diagonal neighbors.

Output File Parameters
----------------------

The **output_name** in the configuration file specifies the name of the output file. The final name depends on the operation and parameters used. You can also specify the output file using the **-output_image** command-line option (e.g., `-output_image toto_out.png`).

For HDF5 output files, the dataset name defaults to `disccofan`.

Pre-processing Parameters
-------------------------

You can set the **preprocessing** option in the config file or via command line to preprocess data. Available options include `ubyte` (8-bit scaling), `ushort` (16-bit scaling), `uint` (32-bit scaling), `log` (natural logarithm), `log10` (log base 10), `exp` (exponential), and `sqrt` (square root).

Although :bold-smallcaps:`disccofan` primarily handles grayscale datasets, it has options for RGB (converted to luminance), LOFAR, or specific sequences. Set **image_options** to "grayscale" (default), "rgb", "lofar", or "sequence" if applicable.

Operation parameters
--------------------

:bold-smallcaps:`disccofan` is one of the most important parameter. It specifies which operation you want to perform on your data. There are 8 possibile operations for analyzing your data: `filter`, `extract`, `dmp`, `csl`, `pattern`, `pattern2D`, `tree`, and `check`, which are described below. Each operation needs some parameters, especially, an attribute function which will specify the morphological properties of the objects you are interested in your data. We first describe what these attributes are, and then the different operations you can use


Attributes
~~~~~~~~~~

There exist 18 attribute functions categorized into 5 groups. Each function is indexed for easy reference and can be used to analyze the morphological properties of objects in the image.


1. **Group #1: Connected Components Size**

   The first group of attributes only contains a single attribute function which determine the size of the objects in your image. It is the simplest property you can analyze with :bold-smallcaps:`disccofan`. The corresponding index of this attribute function is 0.

2. **Group #2: Connected Components Dimensions**

   This group of attributes provides detailed information about the dimensions of objects within your image. Each attribute function returns specific measurements related to the spatial extent of the objects.

   - **Area of Minimal Enclosing Rectangle:** Returns the area of the smallest rectangle that completely encloses the object  (Index 1).
   - **Square of Diagonal of Minimal Enclosing Rectangle:** Provides the square of the diagonal length of the minimal enclosing rectangle (Index 2).
   - **X-Extent:** Returns the maximum extent of the object along the x-axis (Index 3).
   - **Y-Extent:** Returns the maximum extent of the object along the y-axis (Index 4).
   - **Z-Extent:** Returns the maximum extent of the object along the z-axis (Index 5).

   \   
3. **Group #3: Connected Components Center of Mass**

   This group of attributes provides detailed information about the center of mass of the objects within your image. 

   - **X-Mean:** Returns the x coordinates of the center of mass of the object (Index 6).
   - **Y-Mean:** Returns the y coordinates of the center of mass of the object (Index 7).
   - **Z-Mean:** Returns the z coordinates of the center of mass of the object (Index 8).

   \
4. **Group #4: Connected Components Weighted Center of Mass**

   This group of attributes provides detailed information about the center of mass, weighted by the object flux (some of pixel values), of the objects within your image. 

   - **X-Wmean:** Returns the x coordinates of the weighted center of mass of the object (Index 9).
   - **Y-Wmean:** Returns the y coordinates of the weighted center of mass of the object (Index 10).
   - **Z-Wmean:** Returns the z coordinates of the weighted center of mass of the object (Index 11).
   - **Total Flux:** Returns the total flux of the object (Index 12).

   \
5. **Group #5: Connected Components Inertia Properties**

   This group of attributes provides detailed information about the moment of the inertia matrice of each object within your image. 

   - **Trace of the Inertia matrice:** Returns the trace of the inertia matrice correspond to the object (Index 12).
   - **Elongation:** Returns the elongation, which is the ratio of the longest axis over the second longest axis of your object (Index 14).
   - **Flatness:** Returns the flatness (ONLY in 3D!!) which is the ratio of the second longest axis over the thirs longest axis of your object (Index 15).
   - **Sparseness:** Returns the sparseness (also porosity) of the object (Index 16).
   - **Non-compactness:** Returns the non-compactness (also asymettry) of the object (Index 17).

   
Now we will describe the different operations you can do based on these attributes.

1. **Filtering and Extracting Features**: 
   - **`filter`** or **`extract`** apply morphological operations based on chosen attributes. Example syntax:

      .. code-block:: bash
         
         operation = (
            ("filter", 0, 100)
         )
          
   This filters (or extracts if set to extract instead) objects based on their size, with `0` representing the attribute function index and `100` the threshold.

2. **Differential Morphological Profile (DMP)**: 
   - **`dmp`** generates a profile based on a sequence of extract operations with varying thresholds. Example syntax: 

      .. code-block:: bash
   
         operation = (
            ("dmp", 0, "lvec.txt")
         )

   0 sets the attribute function. The file `lvec.txt` should contain the thresholds for filtering.

3. **CSL Segmentation**: 
   - **`csl`** segments the image based on scale, contrast, and luminance. Example syntax: 

      .. code-block:: bash
   
         operation = (
            ("csl", 0, "lvec.txt")
         )

   0 sets the attribute function. The file `lvec.txt` contains thresholds for CSL properties. 
   Note: This operation is partially implemented.

4. **Pattern Spectrum**: 
   - **`pattern`** calculates the 1D pattern spectrum. Example syntax:

      .. code-block:: bash
   
         operation = (
            ("pattern", 0, "lvec.txt", 1)
         )

   0 sets the attribute function. The file `lvec.txt` contains thresholds for defining scales. The scale factor `1` multiplies all values in `lvec.txt`.

5. **2D Pattern Spectrum**: 
   - **`pattern2D`** calculates the 2D pattern spectrum. Example syntax: 

      .. code-block:: bash
   
         operation = (
            ("pattern2D", 0, "lvec.txt", 1, 17, "lvec2.txt", 1)
         )

   The file `lvec.txt` contains the thresholds for the attribute function 0, and `lvec2.txt` contain thresholds for the attribute function 17. The two `1` are the scales used to multiply values in each file. 

6. **Tree**: 
   - **`tree`** generates a text file with the tree structure and attributes. Example syntax: 

      .. code-block:: bash
   
         operation = (
            ("tree", 0)
         )

   Using the tree operation with MPIs will not result in the same output as using it in a sequential manner. Although both approaches are correct, the analysis of a tree file obtained using MPI parallelization must be done carefully, therefore I do not recommend it unless you understand you know exactly what you are doing.

7. **Check**: 
   - **`check`** replaces each pixel value with the corresponding attribute value. Example syntax: 

      .. code-block:: bash
   
         operation = (
            ("check", 0)
         )

   This will fill each pixel with the size (attribute function 0) of the corresponding object it belongs to

Combining several operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to combine several operations at once to minimize the cost of re-building the tree structure and attributes every time. This can be done by setting the operation parameter with a list of several operations, for example:

   .. code-block:: bash
      
      operation = (
         ("filter", 0, 100),
         ("filter", 13, 5),
         ("extract", 5, 1000),
         ("dmp", 1, "tests/lvec.txt"),
         ("pattern", 1, "tests/lvec.txt"),
         ("tree", 1)
      );

All the operations will be performed in a row. If the attribute functions chosen are part of the same group, the computation will be faster since the same tree will be used several times.


Parallelization Options
-----------------------

You can speed up :bold-smallcaps:`disccofan` using threads or MPI. Threads are useful for datasets with a narrow dynamic range (up to 16 bits). Set the number of threads using **threads** (or **-threads**).

For MPI parallelization, use `mpirun` or `mpiexec`. Specify data division with the **mpi_grid** option (or **-grid**). For example, **mpi_grid** set to `[2,2,1]` means 4 MPI processes. The total number of MPI processes must match the grid size.



Tree structure
--------------

By default, :bold-smallcaps:`disccofan` constructs a max-tree, meaning the emphasis is put on the bright objects in the image. If you are interested in the dark objects instead (with the smallest gray values), you can set the option **tree_type** to `min`` for :bold-smallcaps:`disccofan` to build a min-tree.


Verbosity
---------

Set logging level. Options include `off`, `debug`, `info`, `warn`, and `error`.


For help, consult the examples to get a better idea on how to combine these option.



