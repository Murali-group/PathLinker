# Dockerized PathLinker

- Use `COPY` statements to make container more dynamic

To build image from PathLinker repository __root__ directory:

`docker build -t pathlinker/pathlinker -f docker/Dockerfile .`

To run the image interactively:

`docker run -it pathlinker/pathlinker-v3 bash`

To run in Git for Windows:

`winpty docker run -it pathlinker/pathlinker-v3 bash`

To run in Git for Windows and mount the example data inside the container:

`winpty docker run -v /$(pwd)/example:/home/PathLinker/example -it pathlinker/pathlinker-v3 bash`

then inside the container test PathLinker on the example data:

`python run.py --PageRank example/sample-in-net-noweights.txt example/sample-in-nodetypes.txt`
