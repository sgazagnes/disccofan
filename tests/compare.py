import sys

def compare_files(file1, file2):
    """
    Compares two text files line by line and reports the first line where they differ.
    If the files are identical, it prints 'Files are identical'.
    :param file1: Path to the first file.
    :param file2: Path to the second file.
    :return: 0 if files are identical, 1 if they differ
    """
    try:
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            for line_num, (line1, line2) in enumerate(zip(f1, f2), start=1):
                if line1 != line2:
                    print(f"Files differ at line {line_num}")
                    return 1

            # Check if one file has extra lines
            extra_line_f1 = f1.readline()
            extra_line_f2 = f2.readline()

            if extra_line_f1 or extra_line_f2:
                print(f"Files differ at line {line_num + 1}")
                return 1

        print("Files are identical")
        return 0

    except Exception as e:
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare_files.py <file1> <file2>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    sys.exit(compare_files(file1, file2))
