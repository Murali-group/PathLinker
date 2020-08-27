# Dockerized PathLinker

## Steps
- Build PathLinker docker container
- Run PathLinker docker image interactivley and mount example data to container
- Run PathLinker software on example data within the container

## Requirements
- Docker software
- Mac, linux, or Windows machine 

## Build (Mac/Linux/Windows)

To build image: start from PathLinker repository __root__ directory and call:

`docker build -t pathlinker/pathlinker -f docker/Dockerfile .`

To run the image interactively (with no linked data):

`docker run -it pathlinker/pathlinker bash`

## Run (Mac/Linux)

To mount the example data inside the container:

`docker run -v /$(pwd)/example:/home/PathLinker/example -it pathlinker/pathlinker bash`

then inside the container in the `/home/Pathlinker/example` directory call:

`python ../run.py sample-in-net.txt sample-in-nodetypes.txt`

the output files will be linked to the host machine in the example subdirectory. 

To run PathLinker with `PageRank` inside the container:

`python run.py --PageRank example/sample-in-net-noweights.txt example/sample-in-nodetypes.txt`

from the `/home/Pathlinker` directory




## Run (Windows)

To use the docker container with windows, use a bash terminal such as Git for Windows.

To run container in Git for Windows:

`winpty docker run -it pathlinker/pathlinker bash`

**Note the** `winpty` **command**

To run in Git for Windows and mount the example data inside the container:

`winpty docker run -v /$(pwd)/example:/home/PathLinker/example -it pathlinker/pathlinker bash`

then inside the container in the `/home/Pathlinker/example` directory call:

`python ../run.py sample-in-net.txt sample-in-nodetypes.txt`



