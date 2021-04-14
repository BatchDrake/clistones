# clistones
A CLI application to register meteor echoes produced by the french radar system GRAVES

## Requirements
* sigutils (http://github.com/BatchDrake/sigutils) **Note: you will need the most recent version of the develop branch**
* alsa (you need to install libasound2-dev in Debian-based systems)
* CMake 3.11 or newer (although it may work with older versions)

## Building the program
Clone the repo and `cd` to the main directory. Then run:

```
% mkdir build
% cd build
% cmake ..
% make
```
Optionally, you may run `sudo make install` to install it system-wide.

### Help! I'm getting thousand of build errors!
If after running `make` you see errors like these:

```
(metalloid) % make 
Scanning dependencies of target clistones 
[ 33%] Building C object CMakeFiles/clistones.dir/src/graves.c.o 
[ 66%] Building C object CMakeFiles/clistones.dir/src/main.c.o 
In file included from /home/waldo/Documents/Desarrollo/clistones/include/graves.h:28,
                 from /home/waldo/Documents/Desarrollo/clistones/src/graves.c:29:
/usr/local/include/sigutils/sigutils/log.h:74:4: error: #error __FILENAME__ not defined. Please verify your build system.
   74 | #  error __FILENAME__ not defined. Please verify your build system.
      |    ^~~~~
In file included from /home/waldo/Documents/Desarrollo/clistones/src/graves.c:29:
/home/waldo/Documents/Desarrollo/clistones/include/graves.h:115:3: error: unknown type name ‘grow_buf_t’
  115 |   grow_buf_t chirp;
      |   ^~~~~~~~~~

```

This means that you are compiling against a very old version of sigutils, probably because you have cloned the master branch and built from there. You need to build the **develop** branch, which can be cloned by running:

```
% git clone --branch develop https://github.com/BatchDrake/sigutils.git
```

## Running the program
Just plug your radio to the line-in of your computer, tune it to 143.049 kHz USB 
(or 143.051 kHz LSB, although this will reverse the sign of the Doppler) and then 
run `./clistones` (or `clistones` if you installed it system-wide).  You should see
a text line for every echo detected by the program.

## I don't have a radio (yet), how do I test it?
If you have [PulseAudio](https://es.wikipedia.org/wiki/PulseAudio), simply run 
`clistones` as described in the previous step and run `pavucontrol`. In the _Recording_
tab, identify the entry for `ALSA plug-in [clistones]` and click the button on the
right to set the recording device to `Monitor of internal audio` (or something alike).
Then, go to YouTube and [open this video](https://www.youtube.com/watch?v=6T74lSvIc0Y). Make
sure you are listening to it as well. This setup will simulate an actual capture with echoes
recorded during the Perseids meteor shower of 2016.
