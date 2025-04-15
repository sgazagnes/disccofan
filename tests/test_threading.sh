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
echo "| Batch tests: ${bold}THREADING${normal} |"
echo " -------------------------"
echo ""

available_threads=$(nproc)
echo "Found $available_threads processes/threads available"
if [[ $available_threads -lt 4 ]]; then
    echo "Not enough threads for the tests, we won't run the tests with threading."
    exit 1
fi
echo ""

######### Filtering 

echo "1. Testing ${bold}Filtering operations${normal} with multithreading"
echo ""

all_passed=true
run_test "2D Filtering [Area, Non-compactness]" "./disccofan -config tests/config_filtering.ini -output_image tests/out_lena_t4.png --threads 4"
run_compare "4 Threads - Filtered lena - Area " tests/checkimg/GT_area_200.png tests/out_lena_t4_filter_area_200.000.png || all_passed=false
run_compare "4 Threads - Filtered lena - Ncomp " tests/checkimg/GT_ncomp_0.2.png tests/out_lena_t4_filter_ncomp_0.200.png || all_passed=false
display_test_result $all_passed

rm tests/out_lena*

if $all_passed; then
    display_message "32" "${bold}All filtering tests with multithreading passed"
else
    display_message "31" "${bold}Some tests failed"
fi

echo ""
echo ""

######### Attributes 

echo "2. Testing ${bold}Attribute computation${normal} with multithreading"
echo ""
all_passed=true
run_test "2D 4 Threads Attribute [Area, Non-compactness]" "./disccofan -config tests/config_check.ini -output_image tests/out_lena.fits  -threads 4"
run_compare "4 Threads - Attribute 2D lena - Area " tests/checkimg/GT_2D_check_area.fits tests/out_lena_check_area.fits || all_passed=false
run_compare "4 Threads - Attribute 2D lena - Ncomp " tests/checkimg/GT_2D_check_ncomp.fits  tests/out_lena_check_ncomp.fits || all_passed=false
display_test_result $all_passed

all_passed=true
run_test "3D 4 Threads Attribute [Area, Non-compactness]" "./disccofan -config tests/config_check.ini -input_image tests/checkimg/Check16b.fits -output_image tests/out_Check16b.fits  -threads 4"
run_compare "4 Threads - Attribute random 16bits - Area " tests/checkimg/GT_2D_check_area.fits tests/out_lena_check_area.fits || all_passed=false
run_compare "4 Threads - Attribute random 16bits - Ncomp " tests/checkimg/GT_2D_check_ncomp.fits  tests/out_lena_check_ncomp.fits || all_passed=false
display_test_result $all_passed

rm tests/out_Check16b_check*.fits tests/out_lena_check*.fits

if $all_passed; then
    display_message "32" "${bold}All attributes tests with multithreading passed${normal}"
else
    display_message "31" "${bold}Some tests failed${normal}"
fi

echo ""
echo ""

