# clistones
A CLI application to register meteor echoes produced by the french radar system GRAVES

## Requirements
* sigutils (http://github.com/BatchDrake/sigutils)
* alsa (you need to install libasound2-dev in Debian-based systems)
* CMake 3.11 or newer (although it may work with older versions)

## Build
Clone the repo and `cd` to the main directory. Then run:

```
% mkdir build
% cd build
% cmake
% make
```
Optionally, you may run `sudo make install` to install it system-wide.

## Running the program
Just plug your radio to the line-in of your computer, tune it to 145.049 kHz USB 
(or 145.051 kHz LSB, although this will reverse the sign of the Doppler) and then 
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
