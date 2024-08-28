Installation
============

.. role:: bold-smallcaps

To install :bold-smallcaps:`disccofan`, follow these steps:

1. **Clone the Repository**:

   .. code-block:: bash

      git clone https://github.com/sgazagnes/disccofan.git

2. **Install Dependencies**:

   Ensure you have the required libraries installed. :bold-smallcaps:`disccofan` requires the following libraries:
   
   - FreeImage: For reading most image types. (https://freeimage.sourceforge.io/)
   - CFITSIO: For reading FIT/FITS files. (https://heasarc.gsfc.nasa.gov/fitsio/)
   - HDF5: For reading HDF5 (.h5) files. (https://portal.hdfgroup.org/)
   - LibConfig: For reading the config file. (https://hyperrealm.github.io/libconfig/)
   - OpenMP and OpenMPI: For parallelization options

   Some of these libraries may be installed by default on your system. Some libraries are not and you will have to install them manually. 
   
   On a Debian-based system, the following line should make sure you are installing all the necessary libraries:

   .. code-block:: bash

      sudo apt-get install gcc make openmpi-dev libcfitsio-dev libhdf5-dev libfreeimage-dev libconfig-dev 

   On a Fedora or Redhat-based system, the following line should make sure you are installing all the necessary libraries:

   .. code-block:: bash

      sudo yum install gcc make openmpi-devel cfitsio-devel hdf5-devel freeimage-devel libconfig-devel 

   The code was most recently (August 2024) successfully tested on Fedora 40 with the following libraries:

   - `make-1:4.4.1-6.fc40`
   - `gcc-14.2.1-1.fc40`
   - `openmpi-5.0.2-2.fc40`
   - `openmpi-devel-5.0.2-2.fc40`
   - `libconfig-devel-1.7.3-8.fc40`
   - `hdf5-devel-1.12.1-15.fc40`
   - `cfitsio-devel-4.4.0-2.fc40`
   - `freeimage-devel-3.19.0-0.fc40`\

   \

3. **Build the Project**:

   Navigate to the project directory and use `make` to build:

   .. code-block:: bash

      cd disccofan
      make disccofan

3. **Check that the installation worked properly**:

   :bold-smallcaps:`disccofan` has a few bash scripts for testing whether the installation was successfull and that most operations are supported by your system. It is strongly advised to run these tests. Note that these tests require that Python 3 is installed as well as the following packages: numpy, astropy, h5py (run `pip install numpy astropy h5py` if you do not have these packages installed).

   There are three main scripts for testing the code. First one tests the code without any parallelization:

   .. code-block:: bash

       ./tests/test_sequential.sh 

   Second one tests the code with multithreading:

   .. code-block:: bash

       ./tests/test_threading.sh 

   Last one tests the code with MPI parallelization:

   .. code-block:: bash

       ./tests/test_mpi.sh


4. **What to do if the installation failed**:

   If you did not manage to install :bold-smallcaps:`disccofan`, feel free to contact me so we can work a way around. Alternatively, if you know Docker, there exists a Docker image of disccofan that is relatively easy to use so you do not have to go through installing the code and nessecary libraries. See more in the `Running disccofan with Docker` section.