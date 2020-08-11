# Dockerized PathLinker

## DOCKER_BUILD_V1:
- Create and activate Pathlinker environment on docker run call

To build the image from the DOCKER_BUILD_V1 folder:

`docker build -t pathlinker/pathlinker-v1`

To run the image:

`docker run -it pathlinker/pathlinker-v1 bash`

## DOCKER_BUILD_V2:
- activate environment inside image

To build the image from the DOCKER_BUILD_V2 folder:

`docker build -t pathlinker/pathlinker-v2`

To run the image:

`docker run -it pathlinker/pathlinker-v2`

Once the image is running, move into the `/PathLinkerHome/Pathlinker-project/example` diretory. Then call

`bash pl-example.sh` to run pathlinker on the example data.

## DOCKER_BUILD_V3
- use `Copy` statements to make container more dynamic

To build image from Pathlinker repository __root__ directory:

`docker build -t pathlinker/pathlinker-v3 -f docker/DOCKER_BUILD_V3/Dockerfile .`

To run the image:

`docker run -it pathlinker/pathlinker-v3 bash`
