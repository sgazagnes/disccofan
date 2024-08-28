import numpy as np
from astropy.io import fits
from PIL import Image
import h5py
import argparse
import os

def read_image(file_path):
    """
    Reads an image from a file and returns the pixel data as a numpy array.
    Supports FITS, PNG, JPG, and HDF5 formats.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    if file_path.lower().endswith(('.fits', '.fit')):
        with fits.open(file_path) as hdul:
            data = hdul[0].data
            return np.array(data, dtype=np.float32)  # Convert to float32 for consistency
    elif file_path.lower().endswith(('.png', '.jpg', '.jpeg', '.tif', '.tiff')):
        with Image.open(file_path) as img:
            return np.flipud(np.array(img, dtype=np.float32))  # Convert to float32 for consistency
    elif file_path.lower().endswith('.h5'):
        with h5py.File(file_path, 'r') as f:
            dataset_name = list(f.keys())[0]  # Assuming the first dataset is the image data
            data = f[dataset_name][:]
            return np.array(data, dtype=np.float32)  # Convert to float32 for consistency
    else:
        raise ValueError(f"Unsupported file format: {file_path}")

def compare_images(image1, image2):
    """
    Compares two images to determine if they have the same dimensions and pixel values.
    """
    if image1.shape != image2.shape:
        return False, f"Different dimensions: {image1.shape} vs {image2.shape}"

    if not np.array_equal(image1, image2):
        return False, "Pixel values differ"

    return True, "Images are identical"

def main(file1, file2):
    """
    Main function to compare two images and print the result.
    """
    try:
        image1 = read_image(file1)
        image2 = read_image(file2)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return
    except Exception as e:
        print(f"Error reading images: {e}")
        return

    are_same, message = compare_images(image1, image2)
    print(message)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare two images for dimension and pixel value equality.')
    parser.add_argument('file1', type=str, help='Path to the first image file')
    parser.add_argument('file2', type=str, help='Path to the second image file')
    args = parser.parse_args()

    main(args.file1, args.file2)
