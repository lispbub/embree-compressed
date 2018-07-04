###########
##### Setup
###########

In addition to the standard embree requirements our system depends on the eigen3-library.
Before our demos can be started, embree must be set up:
(For further information, we refer to readme.pdf.)

$ cd build
$ cmake ..
$ make -j9

Our code was tested and measurements were taken on Ubuntu 16.04


###########
##### Start
###########

To start a simple interactive demonstration simply run

$ cd build
$ ./viewer -c bomberman.ecs

in the build directory.


###############
##### Box Types
###############

Different leaf types can be enabled by uncommenting the corresponding type in the bomberman.ecs file.
Possible box types are:

1. voxelized leaves
--compress.box

2. pizza box leaves
--compress.leaf

3. full precision grid 
--compress.grid

The subdivision and compression levels can be adjusted with:

# subdivision level (minimum 2)
--subdLvl 6

# compressed levels (minimum 1, maximum 4)
--compLvl 3
