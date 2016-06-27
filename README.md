# Shor
This project is aimed at simulating Shor's algorithm on a classical computer.
All the source code is bundled into one file: Shor.c.

The project can be compiled using the following command:
gcc -o Shor Shor.c -lm
The program is then run from the command line as well.
On unix-based systems, one can test if the program is working by running the program
without any parameters, as such:
./Shor

The program expects two input parameters.
./Shor <number to factor> <accuracy parameter>
The first argument is the number that the program will try to find a divisor of.
If it is not specified, its default value is set to 15.
The second argument is the accuracy parameter. It can be chosen between 0 and 1.
If it is not specified, its default value is set to 0.2.

NOTE: The program will try to find a divisor a maximum number of 100 times. This can
take a long time when trying to factor large numbers though, so first check if the
number is prime before you try to find factors.
It's probably not a good idea to factor numbers larger than 2047 either, as the
program will need more than 4 GB of RAM.

TODO: The program can be optimized by optimizing the function that calculates the
inverse Fourier transform. Finding a better implementation will probably yield
great performance improvements.
