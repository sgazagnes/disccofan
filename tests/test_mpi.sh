#!/bin/bash

bold=$(tput bold)
normal=$(tput sgr0)

# Function to display messages in color
display_message() {
    local color_code=$1
    local message=$2
    echo -e "\033[${color_code}m${message}\033[0m"
}

display_test_result() {
    local test=$1
    if $all_passed; then
        display_message "32" "Test passed"
    else
        display_message "31" "Test failed"
    fi
}

# Function to run compare_images.py and handle output
run_compare() {
    local test_name=$1
    local image1=$2
    local image2=$3
    local output

    output=$(python tests/compare_images.py "$image1" "$image2" 2>&1)
    if [[ $? -ne 0 ]]; then
        #display_message "31" "Test $test_name failed: $output"
        return 1
    elif [[ "$output" == "Images are identical" ]]; then
        #display_message "32" "Test passed"
        return 0
    else
        #display_message "31" "Test $test_name failed: $output"
        return 1
    fi
}

run_compare_txt() {
    local test_name=$1
    local image1=$2
    local image2=$3
    local output

    output=$(python tests/compare.py "$image1" "$image2" 2>&1)
    if [[ $? -ne 0 ]]; then
        #display_message "31" "Test $test_name failed: $output"
        return 1
    elif [[ "$output" == "Files are identical" ]]; then
        #display_message "32" "Test passed"
        return 0
    else
        #display_message "31" "Test $test_name failed: $output"
        return 1
    fi
}


run_test() {
    local test_name=$1
    local command=$2

    echo "Running test: ${bold}$test_name${normal}"
    start_time=$(date +%s)  # Capture the start time
    
    # Run the test
    test_output=$($command 2>&1 | tr -d '\0')
    if [[ $? -ne 0 ]]; then
        display_message "31" "Error running $test_name: $test_output"
        exit 1
    fi

    # Capture the end time and calculate elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    
    echo "Test completed in $elapsed_time seconds"
}

echo ""
echo " ------------------"
echo "| Batch tests: ${bold}MPI${normal} |"
echo " ------------------"
echo ""

# ######### I/O 

echo "1. Testing ${bold}Input/Output${normal}"
echo ""

all_passed=true

run_test "W/O overlap - 4 PNG to 1 FITS" "mpirun -n 4 --oversubscribe ./disccofan -config tests/config_io_parallel.ini -output_image tests/test.fits -tile_overlap 0"
run_compare "W/O overlap - input_image.png vs test.fits" tests/checkimg/input_image.png tests/test.fits || all_passed=false
display_test_result $all_passed
all_passed=true

run_test "W/O overlap - 4 PNG to 1 HDF5" "mpirun -n 4 --oversubscribe ./disccofan -config tests/config_io_parallel.ini -output_image tests/test.h5 -tile_overlap 0"
run_compare "W/O overlap - input_image.png vs test.h5" tests/checkimg/input_image.png tests/test.h5 || all_passed=false
display_test_result $all_passed
rm tests/test.fits tests/test.h5

all_passed=true
run_test "With overlap - 4 PGM to 1 FITS" "mpirun -n 4 --oversubscribe ./disccofan -config tests/config_io_parallel.ini -input_image tests/checkimg/tileov-T0T.pgm -output_image tests/test.fits"
run_compare "Overlap - input_image.png vs test.fits" tests/checkimg/input_image.png tests/test.fits || all_passed=false
display_test_result $all_passed
all_passed=true

run_test "With overlap -4 PGM to 1 HDF5" "mpirun -n 4 --oversubscribe ./disccofan -config tests/config_io_parallel.ini -input_image tests/checkimg/tileov-T0T.pgm -output_image tests/test.h5"
run_compare "Overlap - input_image.png vs test.h5" tests/checkimg/input_image.png tests/test.h5 || all_passed=false
display_test_result $all_passed

all_passed=true
run_test "With overlap - 1 FITS to 4 PNG" "mpirun -np 4 ./disccofan -config tests/config_io_parallel.ini -input_image tests/test.fits -output_image tests/test-T0T.png "
run_compare "Overlap - input_image.png vs test.fits" tests/checkimg/tile-T0T.png tests/test-T0T.png || all_passed=false
display_test_result $all_passed
all_passed=true

run_test "With overlap - 1 HDF5 to 4 FITS" "mpirun -np 4 ./disccofan -config tests/config_io_parallel.ini -input_image tests/test.h5 -output_image tests/test-T0T.fits"
run_compare "Overlap - input_image.png vs test.h5" tests/checkimg/tile-T0T.png tests/test-T0T.fits || all_passed=false
display_test_result $all_passed
rm tests/test.fits tests/test.h5 tests/test-T*T.fits tests/test-T*T.png


if $all_passed; then
    display_message "32" "${bold}All I/O MPI tests passed${normal}"
else
    display_message "31" "${bold}Some tests failed${normal}"
fi

echo ""
echo ""

######## Filtering 

echo "2. Testing ${bold}Filtering operations${normal} with MPI"
echo ""

all_passed=true
run_test "2D 4 MPI Processes Filtering [Area, Non-compactness]" "mpirun -n 4 --oversubscribe ./disccofan -config tests/config_filtering.ini -grid [2,2,1] -output_image tests/out_lena_m4.fits"
run_compare "MPI - Filtered 2D Lena - Area " tests/checkimg/GT_area_200.png tests/out_lena_m4_filter_area_200.000.fits || all_passed=false
run_compare "MPI - Filtered 2D Lena - Ncomp " tests/checkimg/GT_ncomp_0.2.png tests/out_lena_m4_filter_ncomp_0.200.fits || all_passed=false
display_test_result $all_passed

rm tests/out_lena*

if $all_passed; then
    display_message "32" "${bold}All filtering tests passed"
else
    display_message "31" "${bold}Some tests failed"
fi

echo ""
echo ""

######### Attributes 

echo "3. Testing ${bold}Attribute computation${normal}"
echo ""
all_passed=true
run_test "2D 8 MPI Processes Attribute [Area, Non-compactness]" "mpirun -n 8 --oversubscribe ./disccofan  -grid [4,2,1] -config tests/config_check.ini -output_image tests/out_lena.fits"
run_compare "MPI - Attribute Lena - Area " tests/checkimg/GT_2D_check_area.fits tests/out_lena_check_area.fits || all_passed=false
run_compare "MPI - Attribute Lena - Ncomp " tests/checkimg/GT_2D_check_ncomp.fits  tests/out_lena_check_ncomp.fits || all_passed=false
display_test_result $all_passed

all_passed=true
run_test "3D 8 MPI Processes Attribute [Area, Non-compactness]" "mpirun -n 8 --oversubscribe ./disccofan -grid [2,2,2] -config tests/config_check.ini -input_image tests/checkimg/Check16b.fits -output_image tests/out_Check16b.fits "
run_compare "MPI - Attribute Check16b - Area " tests/checkimg/GT_3D_check_area.fits tests/out_Check16b_check_area.fits || all_passed=false
run_compare "MPI - Attribute Check16b - Ncomp " tests/checkimg/GT_3D_check_ncomp.fits tests/out_Check16b_check_ncomp.fits || all_passed=false
display_test_result $all_passed

rm tests/out_Check16b_check*.fits tests/out_lena_check*.fits

if $all_passed; then
    display_message "32" "${bold}All attributes tests passed${normal}"
else
    display_message "31" "${bold}Some tests failed${normal}"
fi

