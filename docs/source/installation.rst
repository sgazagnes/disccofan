Installation
============

.. role:: bold-smallcaps

To install :bold-smallcaps:`disccofan`, follow these steps:

1. **Clone the Repository**:

   .. code-block:: bash

      git clone https://github.com/sgazagnes/disccofan.git

2. **Install Dependencies**:

   Ensure that you have the required libraries installed. :bold-smallcaps:`disccofan` requires the following libraries:
   
   - **FreeImage**: For reading most image types. (https://freeimage.sourceforge.io/)
   - **CFITSIO**: For reading FIT/FITS files. (https://heasarc.gsfc.nasa.gov/fitsio/)
   - **HDF5**: For reading HDF5 (.h5) files. (https://portal.hdfgroup.org/)
   - **LibConfig**: For reading the config file. (https://hyperrealm.github.io/libconfig/)
   - **OpenMP** and **OpenMPI**: For parallelization options

   Some of these libraries may already be installed on your system. If not, you will need to install them manually. 

   On a Debian-based system, use the following command to install the necessary libraries:

   .. code-block:: bash

      sudo apt-get install gcc make openmpi-dev libcfitsio-dev libhdf5-dev libfreeimage-dev libconfig-dev 

   On a Fedora or Redhat-based system, use this command:

   .. code-block:: bash

      sudo yum install gcc make openmpi-devel cfitsio-devel hdf5-devel freeimage-devel libconfig-devel 

   The code was most recently tested (August 2024) on Fedora 40 with the following libraries:

   - `make-1:4.4.1-6.fc40`
   - `gcc-14.2.1-1.fc40`
   - `openmpi-5.0.2-2.fc40`
   - `openmpi-devel-5.0.2-2.fc40`
   - `libconfig-devel-1.7.3-8.fc40`
   - `hdf5-devel-1.12.1-15.fc40`
   - `cfitsio-devel-4.4.0-2.fc40`
   - `freeimage-devel-3.19.0-0.fc40`

3. **Build the Project**:

   Navigate to the project directory and build the project using `make`:

   .. code-block:: bash

      cd disccofan
      make disccofan

4. **Check the Installation**:

   :bold-smallcaps:`disccofan` includes several bash scripts to test the installation and verify that most operations are supported by your system. It is recommended to run these tests. Note that these tests require Python 3 and the following packages: `numpy`, `astropy`, and `h5py`. Install them using:

   .. code-block:: bash

      pip install numpy astropy h5py

   There are three main scripts for testing:

   - To test without parallelization:

     .. code-block:: bash

        ./tests/test_sequential.sh 

   - To test with multithreading:

     .. code-block:: bash

        ./tests/test_threading.sh 

   - To test with MPI parallelization:

     .. code-block:: bash

        ./tests/test_mpi.sh

5. **Troubleshooting Installation Issues**:

   If you encounter problems with the installation, feel free to contact me for assistance. Alternatively, a Docker image of :bold-smallcaps:`disccofan` is available, which simplifies usage without needing to install the code and libraries manually. See the `Running disccofan with Docker` section for more information.
