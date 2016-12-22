//author Daniel Gdanski
//Sharjah, UAE, 19.12.2016

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

double randn (double mu, double sigma);
int any(int array[], int number, int sizeOfArray);
void constraints(int pop[], double gain[], double gamma[], double snr, int numberOfLinks, int chromosomes, int links,
                  int (*cv)[]);
void oneMax(int pop[], double gain[], double gamma[], double snr, int numberOfLinks, int chromosomes, int links,
               double (*fitness)[]);
double meanDouble(double tab[], int size);
double meanInt(int tab[], int size);
double stdDev(double array[], int sizeOfArray);
int xor(int a, int b);
int sumOfArrayElements(int array[], int row, int numberOfElementsInRow);
void geneticAlgorithm(int len, int popSize, int maxGens, int sigmaScalingFlag, int firstParents[], int firstKids[],
	int secondParents[], int secondKids[], int SUSFlag, int crossoverType, int visualizationFlag, int verboseFlag,
	int useMaskRepositoriesFlag, int mutationOnlycrossmasks[], int masks[], int reprodIndices[], int eliteFitness,
	int maskReposFactor, double probCrossover, double probMutation, double sigmaScalingCoeff, double sigma,
	double markers[], int isHybrid, int uniformCrossmaskRepos[], int mutmaskRepos[], int pop[], int eliteIndiv[],
	double cumNormFitnessVals[], double finalFitnessVals[], int finalCvVals[],
	int maxIndexOfFitnessVals, int minIndexOfCvVals, int temp[], int parentIndices[], double avgFitnessHist[],
	double maxFitnessHist[], double avgcvHist[], double mincvHist[], int approxPop[], clock_t start, clock_t stop,
	double gain[], double gamma[], double snr, int numberOfLinks, int approxSol, FILE *f);


int main()
{
    //Generating a random wireless scenario: locations of the transmitters are
    //uniformly distributed ina square area of 1000m x 1000m. The corresponding
    //receiver is randomly located within a radius of 50m and outside a radius
    //of 20m.

    int numberOfLinks = 15;
    int powerOfNumberOfLinks = numberOfLinks * numberOfLinks;
    double snr = pow(10, 6);
    int i, j, k, feasibility, tempInt, indexOfUnsortedGainLink[numberOfLinks], s[numberOfLinks], x[numberOfLinks],
        rows, columns, approxSol;
    double tX[numberOfLinks], tY[numberOfLinks], r[numberOfLinks], theta[numberOfLinks], rX[numberOfLinks],
    rY[numberOfLinks], euclid[powerOfNumberOfLinks], referenceDistance, shadow[powerOfNumberOfLinks], maxShadow, c,
    gain[powerOfNumberOfLinks], gamma[numberOfLinks], gainLink[numberOfLinks], tempDouble, sumS;
    clock_t approxTimeStart, approxTimeStop, startTime, stopTime;

    columns = numberOfLinks;
    rows = 1;
    int approxPop[rows * columns];
    int approxFeasibility[rows];

    //open the .txt file
    FILE *f = fopen("ctext.txt", "w");
    if (f == NULL)
    {
        printf("\nError opening file!\n");
        exit(1);
    }

    srand(time(0));
    for(i = 0; i<numberOfLinks; i++){
        tX[i] = ((double)rand()/(double)RAND_MAX)*1000;
        tY[i] = ((double)rand()/(double)RAND_MAX)*1000;
//        printf("\ntx[%d] = %.3f", i, tX[i]);
        r[i] = 20 + ((double)rand()/(double)RAND_MAX)*30;
//        printf("\nr[%d] = %.3f", i, r[i]);
        theta[i] = ((double)rand()/(double)RAND_MAX)*2* 3.141592653589793;
//        printf("\ntheta[%d] = %.3f", i, theta[i]);

        //interval should be open (0,1)
        if(tX[i] == 0 || tX[i] == 1 || tY[i] == 0 || tY[i] == 1 ||
           r[i] == 0 || r[i] == 30 || theta[i] == 0 || tX[i] == 2* 3.141592653589793){
            i = i - 1;
            continue;
        }

        rX[i] = tX[i] + r[i]*cos(theta[i]);
        rY[i] = tY[i] + r[i]*sin(theta[i]);
    }
    //plot(Tx,Ty,' o')
    //plot(Rx,Ry,' s')

    fprintf(f, "Rx = \n");
    for(i = 0; i<numberOfLinks; i++) fprintf(f, "%.3f ", rX[i]);
    fprintf(f, "\nRy = \n");
    for(i = 0; i<numberOfLinks; i++) fprintf(f, "%.3f ", rY[i]);
    fprintf(f, "\neuclid = \n");

    //Computing the euclidean distance from transmitter of link i to receiver of link j.
    k = 0;
    for(i = 0; i < numberOfLinks; i++){
        for(j = 0; j < numberOfLinks; j++){
            //hypot(x,y) = sqrt(x^2 + y^2)
            euclid[k] = hypot(rX[j] - tX[i], rY[j] - tY[i]);
            fprintf(f, "%.6f ", euclid[k]);
            //printf("\neuclid[%d] = %.6f", k, euclid[k]);
            //find the least link distance
            if (((j == 0 && i == 0 && euclid[k] > 0)) || (euclid[k] < referenceDistance && euclid[k] > 0)){
                //referenceDistance = min{all link distances)
                referenceDistance = euclid[k];
            }
            k++;
        }
        fprintf(f, "\n");
    }

    //Generating random path gains: Assuming log-normal fading;
    //G{i,j}=c*Shadow{i,j}*(euclid{i,j}/dref)^{-1.8} where Shadow_{l} is
    //log-normally distributed with 0 dB mean, and 8 dB
    //standard deviation; c=1/max(Shadow_{i,j})
    for(i = 0; i < powerOfNumberOfLinks; i++){
       //generating A_{l}; log-normal R.V. that represents fading/shadowing
       shadow[i] = pow(10, -(randn(0,1)*2.8284/10));    //8 dB log-variance
       //shadow[i] = pow(10, -(randn(0,1)*1.732/10));    //3 dB log-variance
       if ((i == 0) || (shadow[i] > maxShadow)){
            maxShadow = shadow[i];
       }
       //printf("\nShadow[%d] = %.2f", i,shadow[i]);
    }

    //c = 1/maxShadow;
    c = 0.01;

    //Alternative model: path loss = min[c*Shadow{i,j}*(euclid{i,j}/dref)^{-4},1] for fixed c.
    for(i = 0; i < powerOfNumberOfLinks; i++){
        if((c*shadow[i]*pow(euclid[i]/referenceDistance, -4)) < 1){
            // G(i,j)=min(c*Shadow(i,j)*(euclid(i,j)/dref)^(-4),1);
            gain[i] = (c*shadow[i]*pow(euclid[i]/referenceDistance, -4));
        }else{
            gain[i] = 1;
        }
        //printf("\ngain[%d] = %.9f", i, gain[i]);
    }
    fprintf(f, "G = \n");
    for(i = 0; i < numberOfLinks; i++){
        for(j = 0; j < numberOfLinks; j++){
            fprintf(f, "%.6f ", gain[i*numberOfLinks + j]);
        }
            fprintf(f, "\n");
    }


    //Generate link SINR thresholds at random. Gamma(i) is uniformly distributed
    //over [10 dB, 20 dB]
    for(i = 0; i < numberOfLinks; i++){
        gamma[i] = 10 + ((double)rand()/(double)RAND_MAX)*90;    //SINR threshold ~ U(20 dB, 40 dB)
        //gamma[i] = 100;  //SINR threshold = 20 dB
    }

    //checking problem feasibility
    feasibility = 0;
    for(i = 0; i < numberOfLinks; i++){
        if((snr*gain[i * numberOfLinks + i]) < gamma[i]){
            printf("Problem infeasible");
            feasibility = 1;
            break;
        }
    }

    if (feasibility == 0){
        //Approximation algorithm Implementation

        //tic
        approxTimeStart = clock();

        for(i = 0; i < numberOfLinks; i++){
            gainLink[i] = gain[i * numberOfLinks + i];
            indexOfUnsortedGainLink[i] = i;
            //printf("\nGainlink[%d] [%d] = %.10f", i, i * numberOfLinks + i, gainLink[i]);
        }

        //sort descend gainLink
        for( i = 0; i < numberOfLinks; i++ ){
            for( j = 0; j < numberOfLinks - 1 - i; j++ ){
                if( gainLink[ j ] < gainLink[ j + 1 ] ){
                    tempDouble = gainLink[ j + 1 ];
                    tempInt = indexOfUnsortedGainLink[j+1];
                    gainLink[ j + 1 ] = gainLink[ j ];
                    indexOfUnsortedGainLink[j+1] = indexOfUnsortedGainLink[j];
                    gainLink[ j ] = tempDouble;
                    indexOfUnsortedGainLink[j] = tempInt;
                }
            }
        }
        fprintf(f, "B = \n");
        for( i = 0; i < numberOfLinks; i++ ) fprintf(f, "%.3f ", gainLink[i]);


        for(i = 0; i < numberOfLinks; i++){
            s[i] = -1;
            x[i] = -1;
        }

        fprintf(f, "\nS = \n");
        k = 0;
        for(i = 0; i < numberOfLinks; i++){
            sumS = 0;
            for(j = 0; j < numberOfLinks; j++){
                if((any(s, j, numberOfLinks)) && (j != i) ){
                    //summ_S=summ_S+gamma(I(l))*G(k,I(l))/(G(I(l),I(l))-gamma(I(l))/SNR)+gamma(k)*G(I(l),k)/(G(k,k)-gamma(k)/SNR);
                    sumS = sumS + gamma[indexOfUnsortedGainLink[i]] * (gain[j * numberOfLinks + indexOfUnsortedGainLink[i]]/
                    (gain[indexOfUnsortedGainLink[i]*numberOfLinks + indexOfUnsortedGainLink[i]] - gamma[indexOfUnsortedGainLink[i]]/snr)) +
                    gamma[j] * (gain[indexOfUnsortedGainLink[i] * numberOfLinks + j]/ gain[j*numberOfLinks + j] - gamma[j]/snr);
                }
            }
            if(sumS < 0.5){
                s[k] = indexOfUnsortedGainLink[i];
                fprintf(f, "%d ", s[k]);
                k++;
            }
        }

        fprintf(f, "\nX = \n");
        tempInt = k;
        k = 0;
        for(i = 0; i < tempInt; i++){
            sumS = 0;
            for(j = 0; j < numberOfLinks; j++){
                if((any(s, j, numberOfLinks)) && (j != s[i]) ){
                    //sum2=sum2+gamma(S(i))*G(j,S(i))/(G(S(i),S(i))-gamma(S(i))/SNR);
                    sumS = sumS + gamma[s[i]] * (gain[j * numberOfLinks + s[i]]/ (gain[s[i]*numberOfLinks + s[i]] - gamma[s[i]]/snr));
                }
            }
            if(sumS <= 1){
                x[k] = s[i];
                fprintf(f, "%d ", x[k]);
                k++;
            }
        }


        for(i = 0; i < numberOfLinks; i++){
            approxPop[i] = 0;
            if (any(x, i, numberOfLinks)){
                approxPop[i] = 1;
            }
        }

        fprintf(f, "\napprox_pop = \n");
        for(i = 0; i < columns * rows; i++) fprintf(f, "%d ", approxPop[i]);
        fprintf(f, "\n");

        approxSol = k;
        printf("Approx sol = %d", approxSol);

        printf("\nApprox feasibility = ");
        constraints(approxPop, gain, gamma, snr, numberOfLinks, rows, columns, &approxFeasibility);
        for ( i = 0; i < rows; i++ ) {
            printf( "%d ", approxFeasibility[i]);
        }

        //toc
        approxTimeStop = clock();
        printf("\nApprox time= %.6f", (double) (approxTimeStop - approxTimeStart)/CLOCKS_PER_SEC);


        //Applying (un-hybridized) GA to solve link scheduling
        //
        //SpeedyGA is a vectorized implementation of a Simple Genetic Algorithm in Matlab
        //Version 1.3
        //Copyright (C) 2007, 2008, 2009  Keki Burjorjee
        //Created and tested under Matlab 7 (R14).
        //
        //Licensed under the Apache License, Version 2.0 (the "License"); you may
        //not use this file except in compliance with the License. You may obtain
        //a copy of the License at
        //
        //http://www.apache.org/licenses/LICENSE-2.0
        //
        //Unless required by applicable law or agreed to in writing, software
        //distributed under the License is distributed on an "AS IS" BASIS,
        //WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
        //See the License for the specific language governing permissions and
        //limitations under the License.
        //
        //Acknowledgement of the author (Keki Burjorjee) is requested, but not required,
        //in any publication that presents results obtained by using this script
        //
        //Without Sigma Scaling, Stochastic Universal Sampling, and the generation of mask
        //repositories, SpeedyGA faithfully implements the specification of a simple genetic
        //algorithm given on pages 10,11 of M. Mitchell's book An Introduction to
        //Genetic Algorithms, MIT Press, 1996). Selection is fitness
        //proportionate.

        int len = numberOfLinks;   //The length of the genomes
        int popSize = 50;   //The size of the population (must be an even number)
        int maxGens = 20;   //The maximum number of generations allowed in a run
        double probCrossover = 1;   //The probability of crossing over
        if (probCrossover < 0 || probCrossover > 1){
            printf("\n\nProbability of crossing over is out of range <0,1>\nPut another value.\n");
            return 0;
        }
        double probMutation = 0.003;    //The mutation probability (per bit)
        int sigmaScalingFlag = 1;   //Sigma Scaling is described on pg 168 of M. Mitchell's
                                    //GA book. It often improves GA performance.
        double sigmaScalingCoeff = 1;   //Higher values => less fitness pressure
        double sigma;
        double markers[popSize];
        int firstParents[(popSize/2)*len];
        int firstKids[(popSize/2)*len];
        int secondParents[(popSize/2)*len];
        int secondKids[(popSize/2)*len];
        int SUSFlag = 1;    //1 => Use Stochastic Universal Sampling (pg 168 of M. Mitchell's GA book)
                            //0 => Do not use Stochastic Universal Sampling
                            //Stochastic Universal Sampling almost always improves performance
        int crossoverType = 2;  //0 => no crossover
                                //1 => 1pt crossover
                                //2 => uniform crossover
        int visualizationFlag = 0;  //0 => don't visualize bit frequencies
                                    //1 => visualize bit frequencies
        int verboseFlag = 0;    //1 => display details of each generation
                                //0 => run quietly
        int useMaskRepositoriesFlag = 1;    //1 => draw uniform crossover and mutation masks from
                                            //a pregenerated repository of randomly generated bits.
                                            //Significantly improves the speed of the code with
                                            //no apparent changes in the behavior of
                                            //the SGA
                                            //0 => generate uniform crossover and mutation
                                            //masks on the fly. Slower.

        int mutationOnlycrossmasks[popSize*len];
        int masks[popSize*len];
        int reprodIndices[popSize/2];
        int eliteFitness = -2147483647; //-realmax
        int maskReposFactor = 5;
        int uniformCrossmaskRepos[(popSize/2)*(len + 1)*maskReposFactor];
        int mutmaskRepos[popSize*(len+1)*maskReposFactor];
        int pop[popSize * len];
        int eliteIndiv[len];
        //double fitnessVals[popSize];
        double cumNormFitnessVals[popSize];
        double finalFitnessVals[popSize];
        //int cvVals[popSize];
        int finalCvVals[popSize];
        int maxIndexOfFitnessVals;
        int minIndexOfCvVals = -1;
        int temp[popSize + 1];
        int parentIndices[popSize];

        //preallocate vectors for recording the average and maximum fitness in each generation
        double avgFitnessHist[1*(maxGens+1)], maxFitnessHist[1*(maxGens+1)], avgcvHist[1*(maxGens+1)],
            mincvHist[1*(maxGens+1)];

        int isHybrid = 1;
        geneticAlgorithm(len, popSize, maxGens, sigmaScalingFlag, firstParents, firstKids, secondParents,
                         secondKids, SUSFlag, crossoverType, visualizationFlag, verboseFlag, useMaskRepositoriesFlag,
                         mutationOnlycrossmasks, masks, reprodIndices, eliteFitness, maskReposFactor, probCrossover,
                         probMutation, sigmaScalingCoeff, sigma, markers, isHybrid, uniformCrossmaskRepos, mutmaskRepos,
                         pop, eliteIndiv, cumNormFitnessVals, finalFitnessVals, finalCvVals,
                         maxIndexOfFitnessVals, minIndexOfCvVals, temp, parentIndices, avgFitnessHist, maxFitnessHist,
                         avgcvHist, mincvHist, approxPop, startTime, stopTime, gain, gamma, snr, numberOfLinks, approxSol, f);
        printf("\n");
        isHybrid = 0;
        geneticAlgorithm(len, popSize, maxGens, sigmaScalingFlag, firstParents, firstKids, secondParents,
                         secondKids, SUSFlag, crossoverType, visualizationFlag, verboseFlag, useMaskRepositoriesFlag,
                         mutationOnlycrossmasks, masks, reprodIndices, eliteFitness, maskReposFactor, probCrossover,
                         probMutation, sigmaScalingCoeff, sigma, markers, isHybrid, uniformCrossmaskRepos, mutmaskRepos,
                         pop, eliteIndiv, cumNormFitnessVals, finalFitnessVals, finalCvVals,
                         maxIndexOfFitnessVals, minIndexOfCvVals, temp, parentIndices, avgFitnessHist, maxFitnessHist,
                         avgcvHist, mincvHist, approxPop, startTime, stopTime, gain, gamma, snr, numberOfLinks, approxSol, f);

    }

    //close the .txt file
    fclose(f);

    return 0;
}


void geneticAlgorithm(int len, int popSize, int maxGens, int sigmaScalingFlag, int firstParents[], int firstKids[],
	int secondParents[], int secondKids[], int SUSFlag, int crossoverType, int visualizationFlag, int verboseFlag,
	int useMaskRepositoriesFlag, int mutationOnlycrossmasks[], int masks[], int reprodIndices[], int eliteFitness,
	int maskReposFactor, double probCrossover, double probMutation, double sigmaScalingCoeff, double sigma,
	double markers[], int isHybrid, int uniformCrossmaskRepos[], int mutmaskRepos[], int pop[], int eliteIndiv[],
    double cumNormFitnessVals[], double finalFitnessVals[], int finalCvVals[], int maxIndexOfFitnessVals,
    int minIndexOfCvVals, int temp[], int parentIndices[], double avgFitnessHist[], double maxFitnessHist[],
    double avgcvHist[], double mincvHist[], int approxPop[], clock_t startTime, clock_t stopTime, double gain[],
    double gamma[], double snr, int numberOfLinks, int approxSol, FILE *f){

        int i, j, k, tempInt, gen;
        double sum = 0, tempDouble;
        double fitnessVals[popSize];
        int cvVals[popSize];

        //crossover masks to use if crossoverType==0
		for(i = 0; i < popSize*len; i++){
            mutationOnlycrossmasks[i] = -1;
        }

        //pre-generate two “repositories” of random binary digits from which the
        //the masks used in mutation and uniform crossover will be picked.
        //maskReposFactor determines the size of these repositories.



        for(i = 0; i < popSize*(len+1)*maskReposFactor; i++){
            if (((double)rand()/(double)RAND_MAX) < probMutation){
                mutmaskRepos[i] = 1;
            }else{
                mutmaskRepos[i] = 0;
            }
            //printf("\nmutmask[%d] %d", i, mutmaskRepos[i]);
        }

        fprintf(f, "mutmaskRepos = \n");
        for(i = 0; i < popSize; i++){
            for(j = 0; j < (len+1)*maskReposFactor; j++){
                fprintf(f, "%d ", mutmaskRepos[i * (len+1)*maskReposFactor + j]);
//                printf("%d ", mutmaskRepos[i * (len+1)*maskReposFactor + j]);
            }
                fprintf(f, "\n");
//                printf("\n");
        }

        for(i = 0; i < (popSize/2)*(len + 1)*maskReposFactor; i++){
            if (((double)rand()/(double)RAND_MAX) < 0.5){
                uniformCrossmaskRepos[i] = 1;
            }else{
                uniformCrossmaskRepos[i] = 0;
            }
            //printf("\nuni[%d] %d", i, uniformCrossmaskRepos[i]);
        }


        fprintf(f, "uniformCrossmaskRepos = \n");
        for(i = 0; i < (popSize/2); i++){
            for(j = 0; j < (len+1)*maskReposFactor; j++){
                fprintf(f, "%d ", uniformCrossmaskRepos[i * (len+1)*maskReposFactor + j]);
//                printf("%d ", mutmaskRepos[i * (len+1)*maskReposFactor + j]);
            }
                fprintf(f, "\n");
//                printf("\n");
        }

        //preallocate vectors for recording the average and maximum fitness in each generation
        for(i = 0; i < maxGens+1; i++){
            avgFitnessHist[i] = -1;
            maxFitnessHist[i] = -1;
            avgcvHist[i] = -1;
            mincvHist[i] = -1;
        }

        //the population is a popSize by len matrix of randomly generated boolean values
        for(i = 0; i < popSize * len; i++){
            if (((double)rand()/(double)RAND_MAX) < 0.5){
                pop[i] = 1;
            }else{
                pop[i] = 0;
            }
            //printf("\npop[%d] %d", i, pop[i]);
        }

        if(isHybrid == 1){
            //seeding population with heuristic solution
            for(i = 0; i < popSize/2; i++){
                for(j = 0; j < len; j++){
                    pop[i*len + j] = approxPop[j];
                }
            }
            fprintf(f, "pop = \n");
            for(i = 0; i < popSize; i++){
                for(j = 0; j < len; j++){
                    fprintf(f, "%d ", pop[i * len + j]);
                }
                    fprintf(f, "\n");
            }
        }


        //tic
        startTime = clock();

        for (gen = 0; gen < maxGens; gen++){
            //evaluate the fitness of the population. The vector of fitness values
            //returned  must be of dimensions 1 x popSize.
            for(i = 0; i < popSize; i++) fitnessVals[i] = 0;
            oneMax(pop, gain, gamma, snr, numberOfLinks, popSize, len, &fitnessVals);
            constraints(pop, gain, gamma, snr, numberOfLinks, popSize, len, &cvVals);

            fprintf(f, "fitnessVals = \n");
            for(i = 0; i < popSize; i++) fprintf(f, "%.6f ", fitnessVals[i]);

            fprintf(f, "\ncvVals = \n");
            for(i = 0; i < popSize; i++) fprintf(f, "%d ", cvVals[i]);


            //repository to maintain unscaled/unnormalized fitness and cv
            //values, to access those for the last generation
            maxFitnessHist[gen] = finalFitnessVals[0];
            maxIndexOfFitnessVals = 0;

            for(i = 0; i < popSize; i++){
                finalFitnessVals[i] = fitnessVals[i];
                //find max and index of max
                if(maxFitnessHist[gen] < finalFitnessVals[i]){
                    maxFitnessHist[gen] = finalFitnessVals[i];
                    maxIndexOfFitnessVals = i;
                }
                finalCvVals[i] = cvVals[i];
                //find min and index of min
                if(mincvHist[gen] > finalCvVals[i]){
                    mincvHist[gen] = finalCvVals[i];
                    minIndexOfCvVals  = i;
                }
            }


            avgFitnessHist[gen] = meanDouble(fitnessVals, popSize);
            avgcvHist[gen] = meanInt(cvVals, popSize);


            if(eliteFitness < maxFitnessHist[gen]){
                eliteFitness = maxFitnessHist[gen];
                for(i = 0; i < len; i++){
                    eliteIndiv[i] = pop[maxIndexOfFitnessVals * len + i];
                    //printf("%d ", eliteIndiv[i]);
                }
                //printf("\n");
                fprintf(f, "\neliteIndiv = \n");
                for(i = 0; i < len; i++) fprintf(f, "%d ", eliteIndiv[i]);
            }


            //display the generation number, the average Fitness of the population,
            //and the maximum fitness of any individual in the population
            if (verboseFlag){
                printf("\nGen = %d", gen);
                printf("\navgFitness = %.3f", avgFitnessHist[gen]);
                printf("\nmaxFitness = %.3f", maxFitnessHist[gen]);
            }

            //Conditionally perform bit-frequency visualization
            if(visualizationFlag){
                //here should be plot function
                //figure(1)
                //set (gcf, 'color', 'w');
                //hold off
                //bitFreqs=sum(pop)/popSize;
                //plot(1:len,bitFreqs, '.');
                //axis([0 len 0 1]);
                //title(['Generation = ' num2str(gen) ', Average Fitness = ' sprintf('%0.3f', avgFitnessHist(1,gen+1))]);
                //ylabel('Frequency of the Bit 1');
                //xlabel('Locus');
                //drawnow;
            }

            //Conditionally perform sigma scaling
            if(sigmaScalingFlag){
                sigma = stdDev(fitnessVals, popSize);
                if(sigma != 0){
                    //fitnessVals=1+(fitnessVals-mean(fitnessVals))/(sigmaScalingCoeff*sigma); --MATLAB
                    for(i = 0; i < popSize; i++){
                        fitnessVals[i] = 1 + ((fitnessVals[i] - avgFitnessHist[gen])/(sigmaScalingCoeff*sigma));
                        if(fitnessVals[i] <= 0){
                            fitnessVals[i] = 0;
                        }
                    }
                }else{
                    for(i = 0; i < popSize; i++){
                        fitnessVals[i] = 1;
                    }
                }
            }

            //Normalize the fitness values and then create an array with the
            //cumulative normalized fitness values (the last value in this array
            //will be 1)
            sum = 0;
            for(i = 0; i < popSize; i++){
                sum = sum + fitnessVals[i];
            }
            cumNormFitnessVals[0] = fitnessVals[0]/sum;
            for(i = 1; i < popSize; i++){
                cumNormFitnessVals[i] = cumNormFitnessVals[i - 1] + (fitnessVals[i]/sum);
                j++;
            }
            fprintf(f, "\ncumNormFitnessVals = \n");
            for(i = 0; i < popSize; i++) fprintf(f, "%.3f ", cumNormFitnessVals[i]);


            //Use fitness proportional selection with Stochastic Universal or Roulette
            //Wheel Sampling to determine the indices of the parents
            //of all crossover operations
            if(SUSFlag){
                for(i = 0; i < popSize; i++){
                    markers[i] = ((double)rand()/(double)RAND_MAX) + (i/popSize);
                    if(markers[i] > 1){
                        markers[i] = markers[i] - 1;
                    }
                }
            }else{
                for(i = 0; i < popSize; i++){
                    markers[i] = ((double)rand()/(double)RAND_MAX);
                }
            }

            
            //[temp parentIndices]=histc(markers,[0 cumNormFitnessVals]); --MATLAB
            tempInt = 0;
            for(i = 0; i < popSize + 2; i++){
                temp[i] = 0;
                if(i == 0){
                    for(j = 0; j < popSize; j++){
                        if(markers[j] >= 0 && markers[j] < cumNormFitnessVals[i]){
                            temp[i] = temp[i] + 1;
                            parentIndices[j] = i;
                        }
                    }
                }else{
                    for(j = 0; j < popSize; j++){
                        if(markers[j] >= cumNormFitnessVals[i-1] && markers[j] < cumNormFitnessVals[i]){
                            temp[i] = temp[i] + 1;
                            parentIndices[j] = i;
                        }
                    }
                }
            }

            //randperm
            for(i=0; i<popSize; i++) {
                j = rand()%(popSize-i)+i;
                k = parentIndices[j];
                parentIndices[j] = parentIndices[i];
                parentIndices[i] = k;
            }


            //determine the first parents of each mating pair
            for (i = 0; i < (popSize/2); i++) {
                for(j = 0; j < len; j++){
                    firstParents[i*len + j] = pop[(parentIndices[i]*len) + j];
                }
            }

            //determine the second parents of each mating pair
            j = 0;
            for (i = (popSize/2); i < popSize; i++) {
                for(k = 0; k < len; k++){
                    secondParents[j*len + k] = pop[(parentIndices[i]*len) + k];
                }
                j++;
            }

            fprintf(f, "\nfirstParents = \n");
            for(i = 0; i < (popSize/2); i++){
                for(j = 0; j < len; j++){
                    fprintf(f, "%d ", firstParents[i * len + j]);
                }
                fprintf(f, "\n");
            }

            fprintf(f, "secondParents = \n");
            for(i = 0; i < (popSize/2); i++){
                for(j = 0; j < len; j++){
                    fprintf(f, "%d ", secondParents[i * len + j]);
                }
                fprintf(f, "\n");
            }

            //create crossover masks
            if(crossoverType == 0){
                for(i = 1; i < (popSize*len); i++) masks[i] = mutationOnlycrossmasks[i];
            }else if(crossoverType == 1){
                for(i = 0; i < ((popSize/2)*len); i++) masks[i] = -1;
                for(i = 0; i < (popSize/2); i++){
                    tempDouble = ((double)rand()/(double)RAND_MAX)*(len-1);
                    //temp=ceil(rand(popSize/2,1)*(len-1)); --MATLAb
                    tempInt = tempDouble;
                    if(tempInt <= tempDouble){
                        tempDouble = tempDouble + 1;
                    }else{
                        tempDouble = tempInt;
                    }
                    temp[i] = tempDouble;
                }
                for(i = 0; i < (popSize/2); i++){
                    for(j = 0; j < (temp[i]); j++) masks[i*len + j] = 1;
                }
            }else{
                if(useMaskRepositoriesFlag){
                    //temp=floor(rand*len*(maskReposFactor-1)); --MATLAB
                    tempDouble = ((double)rand()/(double)RAND_MAX)*len*(maskReposFactor-1);
                    //floor tempDouble
                    tempInt = tempDouble;
                    if(tempInt > tempDouble){
                        tempDouble = tempDouble - 1;
                    }else{
                        tempDouble = tempInt;
                    }
                    tempInt = tempDouble;
                    //masks=uniformCrossmaskRepos(:,temp+1:temp+len); --MATLAB
                    k = 0;
                    for(i = 0; i < popSize/2; i++){
                        for(j = tempInt; j < tempInt + len; j++){
                            masks[k] = uniformCrossmaskRepos[i * ((len + 1)*maskReposFactor) + j];
                            k++;
                        }
                    }
                }else{
                    //masks=rand(popSize/2, len)<.5; --MATLAB
                    for(i = 0; i < (popSize/2)*len; i++){
                        masks[i] = ((double)rand()/(double)RAND_MAX);
                        if(masks[i] < 0.5){
                            masks[i] = 1;
                        }else{
                            masks[i] = 0;
                        }
                    }
                }
            }

            //determine which parent pairs to leave uncrossed
            //reprodIndices=rand(popSize/2,1)<1-probCrossover; -- MATLAB
            for(i = 0; i < (popSize/2); i++){
                reprodIndices[i] = ((double)rand()/(double)RAND_MAX);
                if(reprodIndices[i] < 1 - probCrossover){
                    reprodIndices[i] = 1;
                }else{
                    reprodIndices[i] = 0;
                }
            }

            //masks(reprodIndices,:)=false; --MATLAB
            for(i = 0; i < (popSize/2); i++){
                for(j = 0; j < len; j++){
                    if(reprodIndices[i] > 0){
                        masks[(reprodIndices[i] - 1)*len + j] = -1;
                    }
                }
            }

            //implement crossover
            //firstKids=firstParents; --MATLAB
            //secondKids=secondParents; --MATLAB
            for(i = 0; i < (popSize/2)*len; i++){
                firstKids[i] = firstParents[i];
                secondKids[i] = secondParents[i];
            }

            //firstKids(masks)=secondParents(masks); --MATLAB
            //secondKids(masks)=firstParents(masks); --MATLAB
            for(i = 0; i < ((popSize/2)*len); i++){
                if(masks[i] == 1){
                    firstKids[i] = secondParents[i];
                    secondKids[i] = firstParents[i];
                }
            }

            //pop=[firstKids; secondKids]; -- MATLAB
            for(i = 0; i < ((popSize/2)*len); i++){
                pop[i] = firstKids[i];
            }
            j = 0;
            for(i = ((popSize/2)*len); i < 2*((popSize/2)*len); i++){
                pop[i] = secondKids[j];
                j++;
            }

            //implement mutation
            if(useMaskRepositoriesFlag){
                tempDouble = ((double)rand()/(double)RAND_MAX)*len*(maskReposFactor-1);
                //temp=floor(rand*len*(maskReposFactor-1)); --MATLAB
                tempInt = tempDouble;
                if(tempInt > tempDouble){
                    tempDouble = tempDouble - 1;
                }else{
                    tempDouble = tempInt;
                }
                tempInt = tempDouble;          
                //masks=mutmaskRepos(:,temp+1:temp+len); --MATLAB
                k = 0;
                for(i = 0; i < popSize; i++){
                    for(j = tempInt; j < tempInt + len; j++){
                        masks[k] = mutmaskRepos[i * ((len + 1)*maskReposFactor) + j];
                        k++;
                    }
                }
            }else{
                //masks=rand(popSize, len)<probMutation; --MATLAB
                for(i = 0; i < popSize*len; i++){
                    masks[i] = ((double)rand()/(double)RAND_MAX);
                    if(masks[i] < probMutation){
                        masks[i] = 1;
                    }else{
                        masks[i] = 0;
                    }
                }                
            }

            for(i = 0; i < popSize*len; i++){
                pop[i] = xor(pop[i], masks[i]);
            }
            fprintf(f, "pop after xor = \n");
            for(i = 0; i < popSize; i++){
                for(j = 0; j < len; j++){
                    fprintf(f, "%d ", pop[i * len + j]);
                }
                    fprintf(f, "\n");
            }
        }

        if (verboseFlag){
            //here should be plot function
            //figure(2)
            //%set(gcf,'Color','w');
            //hold off
            //%plot([0:maxGens],avgFitnessHist,'k-');
            //hold on
            //plot([0:maxGens],maxFitnessHist,'c-');
            //%hold on
            //plot([0:maxGens],mincvHist,'r-');
            //%plot([0:maxGens],avgcvHist,'r-');
            //title('Maximum and Average Fitness')
            //xlabel('Generation')
            //ylabel('Fitness')
        }

        //toc
        stopTime = clock();

        if(isHybrid == 1){
            printf("\n\nGA time hybrid = %.6f", (double)((stopTime - startTime)/ CLOCKS_PER_SEC));
            //printf("\nmax index in first loop = %d", maxIndexOfFitnessVals);
            printf("\nNumber scheduled links hybrid = %d", sumOfArrayElements(pop, maxIndexOfFitnessVals, len));
            printf("\nFeasibility hybrid = %d", finalCvVals[maxIndexOfFitnessVals]);
            printf("\nImprovement hybrid = %.3f", (double) 100 * (sumOfArrayElements(pop, maxIndexOfFitnessVals, len)-approxSol)/approxSol);
        }else{
            printf("\n\nGA time = %.6f", (double)((stopTime - startTime)/ CLOCKS_PER_SEC));
            //printf("\nmax index in second loop = %d", maxIndexOfFitnessVals);
            printf("\nNumber scheduled links = %d", sumOfArrayElements(pop, maxIndexOfFitnessVals, len));
            printf("\nFeasibility = %d", finalCvVals[maxIndexOfFitnessVals]);
            printf("\nImprovement = %.3f\n", (double) 100 * (sumOfArrayElements(pop, maxIndexOfFitnessVals, len)-approxSol)/approxSol);
        }
	}

	//randn is function which generates random numbers from normal distribution
//full description: https://phoxis.org/2013/05/04/generating-random-numbers-from-normal-distribution-in-c/
//mu - mean 0, sigma 1
double randn (double mu, double sigma){
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}

//if any element of array is equal to number, return 1
int any(int array[], int number, int sizeOfArray){
    int i;
    for(i = 0; i < sizeOfArray; i++){
        if (array[i] == number){
            return 1;
        }
    }
    return 0;
}

//This function sums the entries of each chromosome in pop (MS)
void constraints(int pop[], double gain[], double gamma[], double snr, int numberOfLinks, int chromosomes, int links,
                  int (*cv)[]){
    //chromosomes - rows
    //links - columns
    int i, j, k;
    double interference[chromosomes*links];

    for(i = 0; i < chromosomes; i++){
        (*cv)[i] = 0;
        for(j = 0; j < links; j++){
            interference[i * links + j] = 0;
            for(k = 0; k < links; k++){
                if(k != j){
                    interference[i * links + j] = interference[i * links + j] + pop[i * links + k] * gain[k * numberOfLinks + j];
                }
            }
            if (((gain[j * numberOfLinks + j]/(interference[i * links + j] + 1/snr)) < gamma[j]) && (pop[i * links + j] == 1)){
                (*cv)[i] = (*cv)[i] + 1;
            }
        }
    }
}

//This function sums the entries of each chromosome in pop (MS)
void oneMax(int pop[], double gain[], double gamma[], double snr, int numberOfLinks, int chromosomes, int links,
                double (*fitness)[]){
    //chromosomes - rows
    //links - columns
    double interference[chromosomes*links], cv[chromosomes];
    int i, j, k;
    for(i = 0; i < chromosomes; i++){
        (*fitness)[i] = 0;
        cv[i] = 0;
    }

    for(i = 0; i < chromosomes; i++){
        for(j = 0; j < links; j++){
            (*fitness)[i] = (*fitness)[i] + pop[i * links + j];
        }
        //printf("%.0f  ", (*fitness)[i]);
    }

    for(i = 0; i < chromosomes; i++){
        cv[i] = 0;
        for(j = 0; j < links; j++){
            interference[i * links + j] = 0;
            for(k = 0; k < links; k++){
                if(k != j){
                    interference[i * links + j] = interference[i * links + j] + pop[i * links + k] * gain[k * numberOfLinks + j];
                }
            }
            if (((gain[j * numberOfLinks + j]/(interference[i * links + j] + 1/snr)) < gamma[j]) && (pop[i * links + j] == 1)){
                cv[i] = cv[i] + 1;
                //printf("\ncv[%d]= %.3f", i, cv[i]);
            }
        }
    }

    for(i = 0; i < chromosomes; i++){
        (*fitness)[i] =  (*fitness)[i] / (cv[i] + (double) 1/numberOfLinks);
    }

}

//counts mean value in double
double meanDouble(double tab[], int size){
    int i;
    double mean = 0;
    for(i = 0; i < size; i++){
        mean = mean + tab[i];
    }
    return mean = (double) mean/ (double) size;
}

//counts mean value in int
double meanInt(int tab[], int size){
    int i;
    double mean = 0;
    for(i = 0; i < size; i++){
        mean = mean + tab[i];
    }
    return mean = (double) mean/ (double) size;
}

//standard deviation
double stdDev(double array[], int sizeOfArray){
    double sum = 0, mean;
    int i;
    for(i = 0; i < sizeOfArray; i++){
        sum = sum + array[i];
    }
    mean = sum/sizeOfArray;
    sum = 0;
    for(i = 0; i < sizeOfArray; i++){
        sum = sum + pow(array[i] - mean, 2);
    }
    return sqrt(sum/(sizeOfArray-1));
}

//exclusive or
int xor(int a, int b){
    if(a != b){
        return 1;
    }else{
        return 0;
    }
}

//sum of each element in row in array
int sumOfArrayElements(int array[], int row, int numberOfElementsInRow){
    int i, sum = 0;
    for(i = 0; i < numberOfElementsInRow; i++){
        sum = sum + array[row * numberOfElementsInRow + i];
    }
    return sum;
}
