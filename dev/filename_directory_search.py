import re

def filename_directory_search(input_file, directory):
    """Gives a filename without the path and ensures the directory given by the user is in the expected format."""
    if len(re.split(r'^/(.*)/', input_file)) == 1:
        name = re.match(r'^.*(?=.bam)', input_file).group()
    if len(re.split(r'^/(.*)/', input_file)) > 1:
        no_dir = re.split(r'^/(.*)/', input_file)[-1]
        #input_dir = re.match(r'^/(.*)/', input_file).group()
        name = re.match(r'^.*(?=.bam)', no_dir).group()
    if list(directory)[-1]!="/":
        directory = directory + "/"
    return name, directory
