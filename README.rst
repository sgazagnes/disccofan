DISCCOFAN
=========


.. image:: docs/source/images/disccofan.png
   :width: 300px
   :align: center

|

**disccofan** is a C image processing tool that provides utilities for analyzing patterns and structures in 2D and 3D datasets.

Read the documentation here:

https://disccofan.readthedocs.io/

Techniques that **disccofan** is using were published in these papers:

- `Concurrent Computation of Attribute Filters on Shared Memory Parallel Machines <https://ieeexplore.ieee.org/document/4407727>`_ 
- `Distributed Component Forests in 2-D: Hierarchical Image Representations Suitable for Tera-Scale Images <https://www.worldscientific.com/doi/10.1142/S0218001419400123?srsltid=AfmBOorzh_s6u-6cin0VpWfJYVFr3kvkKw8Chr1SxBPhBGWghmRcXMPG>`_ 
- `Distributed Connected Component Filtering and Analysis in 2D and 3D Tera-Scale Data Sets <https://ieeexplore.ieee.org/document/9376636>`_ 
- `Parallel Attribute Computation for Distributed Component Forests <https://ieeexplore.ieee.org/document/9897660>`_ 

Scaling performance:
----------------
The following pictures present the main scaling results, for running **disccofan** with and without attribute computation, 
comparing its performance with the only other available code. The tests are based on 2-D and 3-D remote-sensing and astronomical datasets with sizes from 5 to 160 Gigapixels.


1. 8 bits per pixel - decreasing tile experiment - MPI + threading:

   .. raw:: html

      <div style="display: flex; flex-direction: column; align-items: center;">
         <img src="docs/scaling/8bit_hybrid.png" width="95%" />
         <p style="text-align: center; font-style: italic; margin-top: 0.5em;">
            Figure 1: Hybrid performance of DISCCOFAN (with only parent-child computation or with area attribute) and GÖTZ-2D on the 2D 8-bpp quantization of the ESO luminance channel (≈ 9 Gpx).
            The size of the individual data chunks is decreasing as the number of MPI processes increases. DISCCOFAN has a small memory overhead, but is faster and scales almost linearly
            up to 64 processes.
         </p>
      </div>

2. floating point - 32 bits per pixel - decreasing tile experiment - MPI + threading:

   .. raw:: html

      <div style="display: flex; flex-direction: column; align-items: center;">
         <img src="docs/scaling/floating_hybrid.png" width="95%" />
         <p style="text-align: center; font-style: italic; margin-top: 0.5em;">
            Figure 2: Hybrid performance of DISCCOFAN (with only parent-child computation or with area attribute) and GÖTZ-2D on the single-precision floating point ESO luminance channel (≈ 9 Gpx).
            The size of the individual data chunks is decreasing as the number of MPI processes increases. 
         </p>
      </div>

3. floating point - 32 bits per pixel - decreasing tile experiment - MPI only:

   .. raw:: html

      <div style="display: flex; flex-direction: column; align-items: center;">
         <img src="docs/scaling/floating_mpi.png" width="30%" />
         <p style="text-align: center; font-style: italic; margin-top: 0.5em;">
            Figure 3: Speed-up of DISCCOFAN and GÖTZ-2 D on the 2D, floating point, ESO luminance channel (≈ 9 Gpx) when using only MPI processes without threads.
         </p>
      </div>

4. floating point - 32 bits per pixel - increasing tile experiment - 162 Gvoxels - MPI only:

   .. raw:: html

      <div style="display: flex; flex-direction: column; align-items: center;">
         <img src="docs/scaling/lofar.png" width="30%" />
         <p style="text-align: center; font-style: italic; margin-top: 0.5em;">
            Figure 4: Execution time (left), speed-up (middle) and memory gain (right) on the LOFAR 3D observation, with single precision floating point values. The 3D
            tile size remains constant (15000 × 15000 × 15 pixels), such that the volume size increases linearly with the number of processes used. With 48 processes,
            the total volume processed is 162 Gvoxels.
         </p>
      </div>


AUTHOR
------

- Simon Gazagnes <sgsgazagnes@gmail.com>

LICENSE
-------

This project is licensed under the MIT License - see the LICENSE file for details.

