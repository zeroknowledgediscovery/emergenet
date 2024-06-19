import gzip
import shutil
import os

def unzip_files(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.gz'):
                gz_file_path = os.path.join(root, file)
                output_file_path = os.path.join(root, file[:-3])
                with gzip.open(gz_file_path, 'rb') as f_in:
                    with open(output_file_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(gz_file_path)

if __name__ == '__main__':
    unzip_files(os.path.dirname(__file__))
