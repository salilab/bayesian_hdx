# Docker image creation for v2

## Build

Download the **Dockerfile** in this folder to the docker server and place it in an empty folder.

Build the image using this command line:

    docker build -t bayesian_hdx .

The docker image place the **bayesian_hdx** code in the folder **/usr/local/bayesian_hdx**. The **PYTHONPATH** 
variable is defined in the image to the folder **/usr/local/bayesian_hdx/pyext/src**

The folder **/usr/local/bayesian_hdx/bin** is added to the global **PATH**. 

The image uses user **ubuntu** and the working directory in the container is **/data**. 
        
## Usage for processing HDXWorkbench CSV file

    docker run -v <path-to-the-directory-with-the-csv>:/data bayesian_hdx hdxworkbench.py -w HDXWorkbench.csv -o
     outputdir

For all options in the script:

    docker run bayesian_hdx hdxworkbench.py --help
    
