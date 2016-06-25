//Standard includes
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"

//Enable the use of booleans
#define true 1
#define false 0
typedef char bool;

/* Greatest Common Divisor */
unsigned int gcd(unsigned int a, unsigned int b)
{
	//Euclidean algorithm
	unsigned int c;
	while (a != 0)
	{
		c = a;
		a = b%a;
		b = c;
	}
	return b;
}

/* Calculcate the remainder of a power with positive offset */
unsigned int remmod(unsigned int base, unsigned int exp, unsigned int offset, unsigned int mod)
{
	//Modular exponentiotion
	unsigned int result = 1;
	base %= mod;
	while (exp != 0)
	{
		if (exp & 1)
		{
			result = (result * base) % mod;
		}
		exp >>= 1;
		base = (base * base) % mod;
	}
	result = (result + offset) % mod;
	return result;
}

/* Complex variables */
typedef struct complex {
	double Re;
	double Im;
} COMPLEX;

/* Inverse Fast Fourier Transform */
void ifftstep(COMPLEX * v, unsigned int N, unsigned int i, COMPLEX * tmp, COMPLEX * rootsofunity)
{
	COMPLEX * ve, * vo;
	COMPLEX z,w;
	if (N > 1)
	{
		ve = tmp;
		vo = &tmp[N/2];
		for (unsigned int k = 0; k < N/2; k++)
		{
			ve[k] = v[2*k];
			vo[k] = v[2*k+1];
		}
		ifftstep(ve,N/2,i<<1,v,rootsofunity);
		ifftstep(vo,N/2,i<<1,v,rootsofunity);
		for (unsigned int k = 0; k < N/2; k++)
		{
			w.Re = rootsofunity[k*i].Re;
			w.Im = rootsofunity[k*i].Im;
			z.Re = w.Re * vo[k].Re - w.Im * vo[k].Im;
			z.Im = w.Re * vo[k].Im + w.Im * vo[k].Re;
			v[k].Re = ve[k].Re + z.Re;
			v[k].Im = ve[k].Im + z.Im;
			v[k+N/2].Re = ve[k].Re - z.Re;
			v[k+N/2].Im = ve[k].Im - z.Im;
		}
	}
	return;
}

void ifft(COMPLEX * v, unsigned int n)
{
	unsigned int N = (1 << n);
	
	//Set up a LUT for the cosines and sines
	COMPLEX * rootsofunity = malloc(N * sizeof(COMPLEX));
	for (unsigned int i = 0; i < N; i++)
	{
		rootsofunity[i].Re = cos(2.0*M_PI*i/(double)N);
		rootsofunity[i].Im = sin(2.0*M_PI*i/(double)N);
	}
	
	COMPLEX * tmp = malloc(N * sizeof(COMPLEX));
	ifftstep(v, N, 1, tmp, rootsofunity);
	
	free(rootsofunity);
	free(tmp);
	return;
}

/* The order finding subroutine */
unsigned int findorder(unsigned int N, double epsilon, unsigned int x)
{
	struct timeval time;
	gettimeofday(&time, NULL);
	srand((unsigned) time.tv_usec * x);
	
	//Calculate the number of qubits needed
	unsigned int L = ceil(log2(N)-1e-8);
	unsigned int t = 2*L + 1 + ceil(log2(2.0+1.0/(2.0*epsilon))-1e-8);
	
	//Allocate space for the calculation and fill with ones
	unsigned int * psi = malloc((1 << t) * sizeof(unsigned int));
	for (unsigned int i = 0; i < (1 << t); i++)
	{
		psi[i] = 1;
	}
	
	//Allocate space for the multiplication matrix and fill it
	unsigned int * multorder = malloc((1 << L) * t * sizeof(unsigned int));
	for (unsigned int i = 0; i < (1 << L); i++)
	{
		if (i >= N) {multorder[i] = i;}
		else {multorder[i] = (i * x) % N;}
	}
	
	for (unsigned int i = 1; i < t; i++)
	{
		for (unsigned int j = 0; j < (1<<L); j++)
		{
			multorder[i*(1<<L)+j] = multorder[(i-1)*(1<<L)+multorder[(i-1)*(1<<L)+j]];
		}
	}
	
	//Apply the multiplications
	unsigned int exp,ctr;
	for (unsigned int i = 0; i < (1 << t); i++)
	{
		exp = i;
		ctr = 0;
		while (exp)
		{
			if (exp & 1)
			{
				psi[i] = multorder[ctr*(1<<L)+psi[i]];
			}
			exp >>= 1;
			ctr++;
		}
	}
	
	//Measure the second register and store the resulting first register in a complex vector	
	unsigned int b = psi[rand() % (1 << t)];
	
	unsigned int sum = 0;
	for (unsigned int i = 0; i < (1 << t); i++) {if (psi[i] == b) {sum++;}}
	double normalization = 1.0/sqrt(sum) / sqrt((double)(1 << t));
	
	COMPLEX * phi = malloc((1 << t) * sizeof(COMPLEX));
	for (unsigned int i = 0; i < (1 << t); i++)
	{
		if (psi[i] == b) {phi[i].Re = normalization;}
		else {phi[i].Re = 0;}
		phi[i].Im = 0.0;
	}
	
	//Apply the inverse Fourier transform
	ifft(phi,t);
	
	//Get the probabilities for measurement
	double p = (double)rand() / (double)RAND_MAX;
	p *= p;
	double cumsum = 0;
	unsigned int k;
	for (k = 0; k < (1 << t); k++)
	{
		cumsum += phi[k].Re * phi[k].Re + phi[k].Im * phi[k].Im;
		if (cumsum >= p) {break;}
	}
	
	//Apply the continuous fraction algorithm
	unsigned int r = 0;
	if (k != 0)
	{
		unsigned int n, d, napprox, dapprox, tmp, ctr;
		int i;
		unsigned int * contfracs = malloc(t * sizeof(unsigned int));
		n = k;
		d = (1 << t);
		ctr = 0;
		while (n != 0)
		{
			contfracs[ctr++] = d/n;
			tmp = n;
			n = d%tmp;
			d = tmp;
		
			napprox = 0;
			dapprox = 1;
			for (i = ctr - 1; i >= 0; i--)
			{
				tmp = contfracs[i] * dapprox + napprox;
				napprox = dapprox;
				dapprox = tmp;
			}
			if (dapprox >= N) {break;}
			r = dapprox;
		}
		
		free(contfracs);
	}
	
	//Deallocate all the memory that was used
	free(psi);
	free(multorder);
	free(phi);
	
	return r;
}

/* Shor's algorithm */
//This structure contains the results of Shor's algorithm
typedef struct Shorstats
{
	unsigned int factor;
	unsigned int x;
	unsigned int r;
	unsigned int errorcode;
	bool quantum;
} SHORSTATS;

void Shor(SHORSTATS * results, unsigned int N, double epsilon, bool msg)
{
	//Set up the random number generator seed
	struct timeval t;
	gettimeofday(&t, NULL);
	srand((unsigned) t.tv_usec);
	
	//Set up the default return values
	results->factor = 0;
	results->x = 0;
	results->r = 0;
	results->errorcode = 0;
	results->quantum = false;
	
	//Check if the number to be factored and epsilon are in the proper range
	if (N <= 2) {results->errorcode = 1; if (msg) {printf("The number is too small to be factored.\n");} return;}
	if (epsilon <= 0 || epsilon >= 1) {results->errorcode = 1; if (msg) {printf("The faulttolerance is invalid.\n");} return;}
	
	//Calculate the effective fault-tolerance
	double faulttolerance = 1.0/(2.0*(pow(2.0,ceil(log2(2.0 + 1.0/(2.0*epsilon)) - 1e-8)) - 2.0));
	
	//Display the welcome message
	if (msg)
	{
		printf("Shor's algorithm initiated...\n");
		printf(" - The number to be factored is %d.\n", N);
		printf(" - The fault-tolerance towards the order finding subroutine is %f.\n",faulttolerance);
	}
	
	//Check if the number to be factored is even
	if (N % 2 == 0) {results->factor = 2; if (msg) {printf("The number is even, so this case is trivial.\nA dividing factor of %d is 2.\n", N);} return;}
	
	//Check if N can be written as a^b with a >= 1 and b >= 2
	double log2N = log2(N);
	for (unsigned int b = 2; b <= log2N; b++)
	{
		double a = pow(2,log2N/b);
		if ((0.5 - fabs(a - floor(a + 1e-8) - 0.5)) <= 1e-8)
		{
			results->factor = floor(a + 1e-8);
			if (msg) {printf("The number can be written as a power.\nA dividing factor of %d is %d.\n", N, results->factor);}
			return;
		}
	}
	
	//Randomly choose an x in the range between 2 and N-2.
	results->x = rand() % (N - 3) + 2;
	if (msg) {printf("The program guesses a number to find the order of. This number is %d.\n", results->x);}
	
	//Return gcd(x,N) if it is a non-trivial divisor of N
	unsigned int temp = gcd(N,results->x);
	if (temp > 1)
	{
		results->factor = temp;
		if (msg) {printf("The program got lucky and guessed a number that is not coprime to N.\nA dividing factor of %d is %d.\n", N, results->factor);}
		return;
	}
	
	//Calculate the order of x modulo N
	results->quantum = true;
	
	unsigned int r1 = findorder(N,epsilon,results->x);
	unsigned int r2 = findorder(N,epsilon,results->x);
	if (r1 == 0 && r2 == 0)
	{
		results->errorcode = 2;
		if (msg) {printf("The order finding subroutine failed.\n");}
		return;
	}
	
	if (r1 == 0) {results->r = r2;}
	else if (r2 == 0) {results->r = r1;}
	else {results->r = r1*r2/gcd(r1,r2);}
	
	//Check if the order is correct
	unsigned int y = 1;
	for (unsigned int a = 0; a < results->r; a++)
	{
		y = (y * results->x) % N;
	}
	if (y != 1) {results->errorcode = 3; if (msg) {printf("The order finding subroutine failed.\n");} return;}
	
	if (msg) {printf("The order is found to be %d.\n",results->r);}
	
	//Check if the order is even
	if (results->r % 2 != 0) {results->errorcode = 4; if (msg) {printf("The order is odd, so the algorithm fails.\n");} return;}
	
	//Check if there is any hope in finding a divisor
	if (remmod(results->x,results->r/2,0,N) == N - 1) {results->errorcode = 5; if (msg) {printf("x^(r/2) = -1 mod N, so the algorithm fails.\n");} return;}
	
	//Return one of the factors
	unsigned int candidate = gcd(remmod(results->x,results->r/2,N-1,N),N);
	if (candidate != 1) {results->factor = candidate; if (msg) {printf("A dividing factor of %d is %d.\n", N, results->factor);} return;}
	candidate = gcd(remmod(results->x,results->r/2,1,N),N);
	if (candidate != 1) {results->factor = candidate; if (msg) {printf("A dividing factor of %d is %d.\n", N, results->factor);} return;}
	
	results->errorcode = 6;
	if (msg) {printf("An unexpected error occured. Exiting.\n");}
	return;
}

/* Main program starts here */
int main(int argc, char * argv[])
{
	//Check the input parameters
	unsigned int N;
	double epsilon;
	if (argc == 1)
	{
		N = 15; printf("No number to factor is supplied, so the default value of 15 is used.\n");
		epsilon = 0.2; printf("No accuracy parameter is supplied, so the default value of 0.2 is used.\n");
	}
	else if (argc == 2) {N = atoi(argv[1]); epsilon = 0.2; printf("No accuracy parameter is supplied, so the default value of 0.2 is used.\n");}
	else if (argc == 3) {N = atoi(argv[1]); epsilon = atof(argv[2]);}
	else {printf("Too many arguments supplied. Expected: ./Shor [<Number to factor>] [<Accuracy parameter>].\n"); return 0;}

	//Set up the statistics
	struct timeval start, stop, starttry, stoptry;
	gettimeofday(&start,NULL);
	SHORSTATS * results = malloc(sizeof(SHORSTATS));
	unsigned int numtries = 0;
	unsigned int numquantumtries = 0;
	unsigned int preliminaryresult = 0;
	
	//Run the algorithm a maximum of 100 times
	for (unsigned int i = 0; i < 100; i++)
	{
		numtries++;
		gettimeofday(&starttry,NULL);
		Shor(results,N,epsilon,false);
		gettimeofday(&stoptry,NULL);
		if (results->errorcode == 0) {preliminaryresult = results->factor;}
		if (results->quantum) {numquantumtries++;}
		if (results->quantum && results->errorcode == 0) {break;}
	}
	gettimeofday(&stop,NULL);
	
	//Display the results
	if (results->quantum && results->errorcode == 0)
	{
		//The program succeeded
		printf("The algorithm succesfully found a factor of %d using the order finding subroutine.\n",N);
		printf("The factor that was found was %d.\n",results->factor);
		printf("The number that was guessed was %d.\n",results->x);
		printf("The order found by the quantum routine was %d.\n",results->r);
		printf("The number of tries was %d.\n",numtries);
		printf("The number of times the order finding routine failed was %d.\n",numquantumtries-1);
		printf("The time elapsed during the succesful try was %f seconds.\n",(stoptry.tv_sec - starttry.tv_sec) + (double)(stoptry.tv_usec - starttry.tv_usec) / 1e6);
		printf("The total elapsed time was %f seconds.\n",(stop.tv_sec - start.tv_sec) + (double)(stop.tv_usec - start.tv_usec) / 1e6);
	}
	else
	{
		//The program failed
		printf("The program did not find a dividing factor of %d using the order finding subroutine after 100 tries.\n",N);
		if (preliminaryresult == 0)
		{
			printf("Neither did it find a dividing factor using classical methods. Perhaps the number is prime.\n");
		}
		else
		{
			printf("However, it did find a dividing factor using classical methods.\n");
			printf("The factor that was found was %d.\n",preliminaryresult);
		}
	}

	free(results);
	return 0;
}
