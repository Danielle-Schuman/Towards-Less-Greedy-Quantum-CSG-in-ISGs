import os
import shutil

def move_files_with_substring(source_dir, substring, destination_dir):
    # Create destination directory if it doesn't exist
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    # Walk through the directory tree
    for root, dirs, files in os.walk(source_dir):
        print(root)
        for filename in files:
            # Check if the filename contains the specified substring
            if substring in filename:
                source_file = os.path.join(root, filename)
                destination_file = os.path.join(destination_dir, filename)
                # Move the file to the destination directory
                shutil.move(source_file, destination_file)
                print(f"Moved {source_file} to {destination_file}")

# Define source directory, substring, and destination directory
source_directory = "dwave/parallel"
substring_to_search = "at_most"
destination_directory = f"{source_directory}/old_at_most"

# Call the function to move files
move_files_with_substring(source_directory, substring_to_search, destination_directory)

print("All files containing 'at_most' substring have been moved to", destination_directory)

