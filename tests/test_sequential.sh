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
    test_output=$($command 2>&1)
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
echo " -------------------------"
echo "| Batch tests: ${bold}SEQUENTIAL${normal} |"
echo " -------------------------"
echo ""

# ######### I/O 

echo "1. Testing ${bold}Input/Output${normal}"
echo ""

all_passed=true

run_test "PNG to FITS" "./disccofan -config tests/config_io.ini -output_image tests/test.fits"
run_compare "FITS to HDF5" tests/checkimg/input_image.png tests/test.fits || all_passed=false
display_test_result $all_passed

run_test "FITS to HDF5" "./disccofan -config tests/config_io.ini -input_image tests/test.fits -output_image tests/test.h5"
run_compare "HDF5 to PNG" tests/checkimg/input_image.png tests/test.h5 || all_passed=false
display_test_result $all_passed

run_test "HDF5 to PNG" "./disccofan -config tests/config_io.ini -input_image tests/test.h5 -output_image tests/test.png"
run_compare "PNG to FITS" tests/checkimg/input_image.png tests/test.png || all_passed=false
display_test_result $all_passed

if $all_passed; then
    display_message "32" "${bold}All I/O tests passed${normal}"
else
    display_message "31" "${bold}Some tests failed${normal}"
fi

rm tests/test.png tests/test.fits tests/test.h5

echo ""
echo ""

######### Filtering 

echo "2. Testing ${bold}Filtering operations${normal}"
echo ""

all_passed=true
run_test "2D Filtering [Area, Non-compactness]" "./disccofan -config tests/config_filtering.ini -output_image tests/out_lena.png"
run_compare "Sequential - Filtered lena - Area " tests/checkimg/GT_area_200.png tests/out_lena_filter_area_200.000.png || all_passed=false
run_compare "Sequential - Filtered lena - Non-compactness " tests/checkimg/GT_ncomp_0.2.png tests/out_lena_filter_ncomp_0.200.png || all_passed=false
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
run_test "2D Attribute Check [Area, Non-compactness]" "./disccofan -config tests/config_check.ini -output_image tests/out_lena.fits"
run_compare "Sequential - Attribute 2D lena - Area " tests/checkimg/GT_2D_check_area.fits tests/out_lena_check_area.fits || all_passed=false
run_compare "Sequential - Attribute 2D lena - Ncomp " tests/checkimg/GT_2D_check_ncomp.fits  tests/out_lena_check_ncomp.fits || all_passed=false
display_test_result $all_passed

all_passed=true
run_test "3D Attribute Check [Area, Non-compactness]" "./disccofan -config tests/config_check.ini -input_image tests/checkimg/Check16b.fits -output_image tests/out_Check16b.fits"
run_compare "Sequential - Attribute 3D random 16bits - Area " tests/checkimg/GT_3D_check_area.fits tests/out_Check16b_check_area.fits || all_passed=false
run_compare "Sequential - Attribute 3D random 16bits - Ncomp " tests/checkimg/GT_3D_check_ncomp.fits tests/out_Check16b_check_ncomp.fits || all_passed=false
display_test_result $all_passed

rm tests/out_Check16b_check*.fits tests/out_lena_check*.fits

if $all_passed; then
    display_message "32" "${bold}All attributes tests passed${normal}"
else
    display_message "31" "${bold}Some tests failed${normal}"
fi

echo ""
echo ""

######### Tree extraction 

echo "4. Testing ${bold}Tree extraction${normal}"
echo ""
all_passed=true
run_test "Tree writing" "./disccofan -operation \(\(tree,15\)\) -output_image tests/lena.txt"
run_compare_txt "Sequential - Tree writing" tests/checkimg/GT_tree.txt tests/lena_tree.txt || all_passed=false
display_test_result $all_passed

rm tests/lena_tree.txt

if $all_passed; then
    display_message "32" "${bold}All tree writing tests passed${normal}"
else
    display_message "31" "${bold}Some tests failed${normal}"
fi
