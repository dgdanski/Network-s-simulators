//author: Daniel Gdanski
//Sharjah,  UAE, 25.10.2016

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

double randn (double mu, double sigma);

int main()
{
    //required tolerance
    double epsilon = pow(10, -5);
    //minimum required spectral efficiency (in bits/s/Hz)
    int minSpectralEfficiency = 1;
    double maxAllowedPower = 10;
    double noisePower = pow(10, -5);
    int numberOfNodes = 20;
    //source nodes
    int source = 1;
    //destination nodes
    int destination = 5;
    double x[numberOfNodes], y[numberOfNodes];
    int i, j, k, m, tt, hh, counter, index, powerOfNumberOfNodes = numberOfNodes * numberOfNodes, orderNumber = 1,
        numberOfZeros[powerOfNumberOfNodes/(numberOfNodes+1)], previous[numberOfNodes][numberOfNodes], rPath1[powerOfNumberOfNodes];
    double pMin, pMax, tempPower, maxShadow, c, rWidest, maxR, referenceDistance = 100,  width[powerOfNumberOfNodes],
        length[powerOfNumberOfNodes], linkWeights[powerOfNumberOfNodes], euclid[numberOfNodes * numberOfNodes],
        shadow[powerOfNumberOfNodes], r[numberOfNodes];
    bool links[powerOfNumberOfNodes];
    clock_t begin, end;

    if (numberOfNodes < destination){
        printf("Number of nodes can't be less than destination!\nChange the value and please try again.\n\n");
        return 0;
    }

    srand(time(0));
    for(i = 0; i < numberOfNodes; i++){
        x[i] = rand()%101 + ((double)rand()/(double)RAND_MAX);
        y[i] = rand()%101 + ((double)rand()/(double)RAND_MAX);
//        printf("x[%d] = %f ", i, x[i]);
//        printf("y[%d] = %f \n", i, y[i]);
    }

    k = 0;
    for(i = 0; i < numberOfNodes; i++){
        for(j = 0; j < numberOfNodes; j++){
                //hypot(x,y) = sqrt(x^2 + y^2)
                euclid[k] = hypot(x[j] - x[i],y[j] - y[i]);
                //find the least link distance
                if (((j == 0 && i == 0 && euclid[k] > 0)) || (euclid[k] < referenceDistance && euclid[k] > 0)){
                    //referenceDistance = min{all link distances)
                    referenceDistance = euclid[k];
                }
                k++;
        }
    }
    //display euclid
//    for(m = 0; m < powerOfNumberOfNodes; m++){
//        printf("\neuclid[%d] = %.4f", m, euclid[m]);
//    }

    numberOfZeros[0] = 1;
    j = 1;
    for(i = 1 + (numberOfNodes+1); i <= powerOfNumberOfNodes; i = i + (numberOfNodes+1)){
        if (i == 1) continue;
            numberOfZeros[j] = i;
            //printf("Zeros are in zero[%d] %d\n", j,i);
            j++;
    }

//Generating the random link weights log(1+PG_{l}/N_{0}B)
//assuming log-normal fading
//G_{l}=c*Shadow_{l}*(d_{l}/dref)^{-4}
//where Shadow_{l} is log-normally distributed with 0 dB mean, and 8 dB standard deviation.
//c=1/max(Shadow_{l})

    j = 0;
    orderNumber = 1;
    for(i = 0; i < powerOfNumberOfNodes; i++){
       if(orderNumber == numberOfZeros[j]){
            links[i] = 0;
            j++;
       }else{
            links[i] = 1;
       }
       //printf("Link[%d] = %d\n", i,links[i]);
       orderNumber++;
       //generating A_{l}; log-normal R.V. that represents fading/shadowing
       shadow[i] = links[i]*pow(10, -(randn(0,1)*2.8284/10));
       if ((i == 0) || (shadow[i] > maxShadow)){
            maxShadow = shadow[i];
       }
       //printf("\nShadow[%d] = %.2f", i,shadow[i]);
    }
    //printf("\nMax shadow is %.2f\n", maxShadow);

    //display shadow
//    for(i = 0; i < powerOfNumberOfNodes; i++){
//        printf("\nShadow[%d] = %.2f ", i, shadow[i]);
//    }

    c = 1/maxShadow;

    //generating the link weights
    for(i = 0; i < powerOfNumberOfNodes; i++){
        if (links[i] == 1){
            linkWeights[i] = log2(1+(maxAllowedPower/noisePower)*c*shadow[i]*pow((euclid[i]/referenceDistance), -4));
        }else{
            linkWeights[i] = 0;
        }
    }

    //display link weights
//    for(m = 0; m < powerOfNumberOfNodes; m++){
//        printf("\nlinkweight[%d] = %.2f", m, linkWeights[m]);
//    }

    //initializing Bellman-Ford algorithm for widest path
    counter = 1;
    for(i = 0; i < powerOfNumberOfNodes; i++){
        if (i%numberOfNodes == 0){
            if (counter == source){
                for(j = i; j < i + numberOfNodes; j++){
                        width[j] = INFINITY;
                        length[j] = 0;
                }
            }else{
                width[i] = 0;
                length[i] = INFINITY;
            }
        counter++;
        }
    }

    //fill 'previous' matrix with '0'
    for(i=0; i<numberOfNodes; i++){
        for(j=0; j<numberOfNodes; j++){
            previous[i][j] = 0;
        }
    }

    r[0] = 0;
    //Bellman-Ford Iterations
    for(i = 1; i < numberOfNodes; i++){
        if (i%numberOfNodes != 0){
            for(m = 0; m < numberOfNodes; m++){
                width[m*numberOfNodes + i] = width[m*numberOfNodes + i - 1];
                length[m*numberOfNodes + i] = length[m*numberOfNodes + i - 1];
            }
        }
        for(j = 0; j < numberOfNodes; j++){
            for(k = 0; k < numberOfNodes; k++){
                //finding min(width,linkWeights) and then...
                if ((width[i-1+(j*numberOfNodes)] < linkWeights[k+j*numberOfNodes]) && (links[k+j*numberOfNodes] == 1)){
                    //if (min(width,linkWeights) = width)>another width then...
                    if (width[i-1+(j*numberOfNodes)] > width[i+k*numberOfNodes]){
                        width[i+k*numberOfNodes] = width[i-1+(j*numberOfNodes)];
                        length[i+k*numberOfNodes] = length[i-1+(j*numberOfNodes)] + 1;
                        previous[k][i] = j+1;
                    }
                }else if((linkWeights[k+j*numberOfNodes] < width[i-1+j*numberOfNodes])&& (links[k+j*numberOfNodes] == 1)){
                    //if (min(width,linkWeights) = linkWeights)>another width then...
                    if (linkWeights[k+j*numberOfNodes] > width[i+k*numberOfNodes]){
                        width[i+k*numberOfNodes] = linkWeights[k+j*numberOfNodes];
                        length[i+k*numberOfNodes] = length[i-1+(j*numberOfNodes)] + 1;
                        previous[k][i] = j+1;
                    }
                //if width == linkWeights then...
                }else{
                    //if (width or linkWeights (this is the same value -> watch condition above))>another width then...
                    if ((linkWeights[k+j*numberOfNodes] > width[i+k*numberOfNodes]) && (links[k+j*numberOfNodes] == 1)){
                        width[i+k*numberOfNodes] = linkWeights[k+j*numberOfNodes];
                        length[i+k*numberOfNodes] = length[i-1+(j*numberOfNodes)] + 1;
                        previous[k][i] = j+1;
                    }
                }

            }
        }
        r[i] = width[(destination-1)*numberOfNodes+(i)]/length[(destination-1)*numberOfNodes+(i)];
   }

    //find max(r)
    maxR = r[0];
    for(i = 0; i < numberOfNodes; i++){
        if (r[i] > maxR){
            maxR = r[i];
            index = i;
        }
    }

    rWidest = maxR;
    rPath1[0] = destination;
    tt = destination;
    hh = index;
    i = 1;
    while (source != tt){
        rPath1[i] = previous[tt-1][hh];
        tt = previous[tt-1][hh];
        hh--;
        i++;
    }

    if(rWidest < minSpectralEfficiency){
        printf("Rpath1 = ");
        for(j = 0; j < i; j++){
            printf("%d ", rPath1[j]);
        }
        printf("\nR_widest = %f", rWidest);
        printf("\nP = %f", maxAllowedPower);
        printf("\nProblem is infeasible\n\n");
    }else{
        //tic
        begin = clock();

        //beginning bisection iterations
        pMin = 0;
        pMax = maxAllowedPower;

        while((pMax - pMin) > epsilon){
            //clear width, length
            for(i=0; i<powerOfNumberOfNodes; i++){
                width[i] = 0;
                length[i] = 0;
            }
            //fill previous matrix with '0' == clear previous and clear r
            for(i=0; i<numberOfNodes; i++){
                for(j=0; j<numberOfNodes; j++){
                    previous[i][j] = 0;
                }
                r[i] = 0;
            }
            tempPower = (pMin + pMax)/2;

            //generating the link weights
            for(i = 0; i < powerOfNumberOfNodes; i++){
                if (links[i] == 1){
                    linkWeights[i] = log2(1+(tempPower/noisePower)*c*shadow[i]*pow((euclid[i]/referenceDistance), -4));
                }else{
                    linkWeights[i] = 0;
                }
            }

            //initializing Bellman-Ford algorithm for widest path
            counter = 1;
            for(i = 0; i < powerOfNumberOfNodes; i++){
                if (i%numberOfNodes == 0){
                    if (counter == source){
                        for(j = i; j < i + numberOfNodes; j++){
                            width[j] = INFINITY;
                            length[j] = 0;
                        }
                    }else{
                    width[i] = 0;
                    length[i] = INFINITY;
                    }
                counter++;

                }

            }

            //Bellman-Ford Iterations
            for(i = 1; i < numberOfNodes; i++){
                if (i%numberOfNodes != 0){
                    for(m = 0; m < numberOfNodes; m++){
                        width[m*numberOfNodes + i] = width[m*numberOfNodes + i - 1];
                        length[m*numberOfNodes + i] = length[m*numberOfNodes + i - 1];
                    }
                }

                for(j = 0; j < numberOfNodes; j++){
                    for(k = 0; k < numberOfNodes; k++){
                        //finding min(width,linkWeights) and then...                        //printf("\nLinks[%d] = %d", k+j*numberOfNodes, links[k+j*numberOfNodes]);
                        if ((width[i-1+(j*numberOfNodes)] < linkWeights[k+j*numberOfNodes]) && (links[k+j*numberOfNodes] == 1)){
                            //if (min(width,linkWeights) = width)>another width then...
                            if (width[i-1+(j*numberOfNodes)] > width[i+k*numberOfNodes]){
                                width[i+k*numberOfNodes] = width[i-1+(j*numberOfNodes)];
                                length[i+k*numberOfNodes] = length[i-1+(j*numberOfNodes)] + 1;
                                previous[k][i] = j+1;
                            }
                        }else if((linkWeights[k+j*numberOfNodes] < width[i-1+j*numberOfNodes])&& (links[k+j*numberOfNodes] == 1)){
                            //if (min(width,linkWeights) = linkWeights)>another width then...
                            if (linkWeights[k+j*numberOfNodes] > width[i+k*numberOfNodes]){
                                width[i+k*numberOfNodes] = linkWeights[k+j*numberOfNodes];
                                length[i+k*numberOfNodes] = length[i-1+(j*numberOfNodes)] + 1;
                                previous[k][i] = j+1;
                            }
                        //if width == linkWeights then...
                        }else{
                            //if (width or linkWeights (this is the same value -> watch condition above))>another width then...
                            if ((linkWeights[k+j*numberOfNodes] > width[i+k*numberOfNodes]) && (links[k+j*numberOfNodes] == 1)){
                                width[i+k*numberOfNodes] = linkWeights[k+j*numberOfNodes];
                                length[i+k*numberOfNodes] = length[i-1+(j*numberOfNodes)] + 1;
                                previous[k][i] = j+1;
                            }
                        }

                    }
                }
                r[i] = width[(destination-1)*numberOfNodes+(i)]/length[(destination-1)*numberOfNodes+(i)];
            }

            //find max(r)
            maxR = r[0];
            for(i = 0; i < numberOfNodes; i++){
                if (r[i] > maxR){
                    maxR = r[i];
                    index = i;
                }
            }
            rWidest = maxR;

            if (rWidest > minSpectralEfficiency){
                pMax = tempPower;
            }else{
                pMin = tempPower;
            }
        }

        //toc
        end = clock();

        rPath1[0] = destination;
        tt = destination;
        hh = index;
        i = 1;
        while (source != tt){
            rPath1[i] = previous[tt-1][hh];
            tt = previous[tt-1][hh];
            hh--;
            i++;
        }
        printf("Elapsed time is = %.3f seconds", ((double)(end - begin) / CLOCKS_PER_SEC));
        printf("\nRpath1 = ");
        for(j = 0; j < i; j++){
            printf("%d ", rPath1[j]);
        }
        printf("\nR_widest = %f", rWidest);
        printf("\nP = %e\n\n", tempPower);
    }
    return 0;
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
