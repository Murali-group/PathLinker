# Dockerized PathLinker

- Use `COPY` statements to make container more dynamic

To build image from PathLinker repository __root__ directory:

`docker build -t pathlinker/pathlinker -f docker/Dockerfile .`

To run the image interactively:

`docker run -it pathlinker/pathlinker bash`

To run in Git for Windows:

`winpty docker run -it pathlinker/pathlinker bash`

To run in Git for Windows and mount the example data inside the container:

`winpty docker run -v /$(pwd)/example:/home/PathLinker/example -it pathlinker/pathlinker bash`

then inside the container in the `/home/Pathlinker/example` directory call:

`python ../run.py sample-in-net.txt sample-in-nodetypes.txt`

the output files will be linked to the host machine in the example subdirectory. 

To run Pathlinker inside the container:

`python run.py --PageRank example/sample-in-net-noweights.txt example/sample-in-nodetypes.txt`

from the `/home/Pathlinker` directory
