//author: Daniel Gdanski
//Sharjah,  UAE, 9.11.2016

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

clock_t timeSum(clock_t tab[], int numberOfElements);
int sumInt(int tab[], int numberOfElements);
double sum(double tab[], int numberOfElements);
double min(double x, double y);
double max(double x, double y);
double randn (double mu, double sigma);

int main()
{
    //minimum required spectral efficiency (in bits/s/Hz)
    double gamma=1;
    //infinity
    int inf = 2147483647;

    double noisePower = pow(10, -12);
    int numberOfNodes = 5;
    int source = 1;
    int destination = numberOfNodes;
    int totalExperiments = 10000;
    int i, j, k, m, counter, exper, index[totalExperiments];
    int orderNumber, powerOfNumberOfNodes = numberOfNodes * numberOfNodes,
        numberOfZeros[powerOfNumberOfNodes/(numberOfNodes+1)];
    double x[numberOfNodes], y[numberOfNodes];
    double c, referenceDistance, maxShadow, pMax[numberOfNodes], euclid[powerOfNumberOfNodes], shadow[powerOfNumberOfNodes],
        linkWeights[powerOfNumberOfNodes], width[powerOfNumberOfNodes], length[powerOfNumberOfNodes], pMaxMin[totalExperiments],
        directLinkPower[totalExperiments];
    bool links[powerOfNumberOfNodes];
    clock_t tic, toc, measuredTime[totalExperiments];

    if (source > numberOfNodes){
        printf("Source node can't be bigger than number of nodes!\nChange the value and please try again.\n\n");
        return 0;
    }

    //fill pMaxMin with inf
    for(i=0; i<totalExperiments; i++){
        pMaxMin[i] = inf;
    }
    for (exper = 0; exper < totalExperiments; exper++){

        //Generating the nodes, assume 20 nodes, randomly scattered in a 100*100 area,
        //a link exists if the distance between two nodes is >= 25 m

        x[0] = 0;
        y[0] = 0;
        x[numberOfNodes-1] = 100;
        y[numberOfNodes-1] = 100;

        srand(time(0) - exper);
        for(i = 1; i < numberOfNodes - 1; i++){
            x[i] = rand()%101 + ((double)rand()/(double)RAND_MAX);
            y[i] = rand()%101 + ((double)rand()/(double)RAND_MAX);
            //printf("\nx[%d] = %.12f", i, x[i]);
            //printf("\ny[%d] = %.12f", i, y[i]);
        }

        referenceDistance = hypot(100, 100);
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
        //reference distance for the far field; here we take referenceDistance = 0.1
        referenceDistance = 0.1;


        numberOfZeros[0] = 1;
        j = 1;
        for(i = 1 + (numberOfNodes+1); i <= powerOfNumberOfNodes; i = i + (numberOfNodes+1)){
            if (i == 1) continue;
            numberOfZeros[j] = i;
            j++;
        }


        //Generating the random link weights log(1+PG_{l}/N_{0}B)
        //assuming log-normal fading
        //G_{l}=c*Shadow_{l}*(d_{l}/dref)^{-4}
        //where Shadow_{l} is log-normally distributed with 0 dB mean, and 8 dB standard deviation.
        //c=1/max(Shadow_{l})
        maxShadow = 0;
        j = 0;
        orderNumber = 1;
        for(i = 0; i < powerOfNumberOfNodes; i++){
           if(orderNumber == numberOfZeros[j]){
                links[i] = 0;
                j++;
           }else{
                links[i] = 1;
           }
           orderNumber++;
           //generating A_{l}; log-normal R.V. that represents fading/shadowing
           shadow[i] = links[i]*pow(10, -(randn(0,1)*2.8284/10));
           if ((i == 0) || (shadow[i] > maxShadow)){
                maxShadow = shadow[i];
           }
        }

        c = 1/maxShadow;
        c = 0.01;

        //generating the link weights
        for(i = 0; i < powerOfNumberOfNodes; i++){
            if (links[i] == 1){
                linkWeights[i] = (1/noisePower)*c*shadow[i]*pow(max(referenceDistance, euclid[i]), -4);
                if(linkWeights[i] > (1/noisePower)){
                    linkWeights[i] = (1/noisePower);
                }
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
                            width[j] = inf;
                            length[j] = 0;
                    }
                }else{
                    width[i] = 0;
                    length[i] = inf;
                }
            counter++;
            }
        }

        //clear pMax
        for(i=0; i<numberOfNodes; i++){
            pMax[i] = 0;
        }

        //tic
        tic = clock();

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
                        }
                    }else if((linkWeights[k+j*numberOfNodes] < width[i-1+j*numberOfNodes])&& (links[k+j*numberOfNodes] == 1)){
                        //if (min(width,linkWeights) = linkWeights)>another width then...
                        if (linkWeights[k+j*numberOfNodes] > width[i+k*numberOfNodes]){
                            width[i+k*numberOfNodes] = linkWeights[k+j*numberOfNodes];
                            length[i+k*numberOfNodes] = length[i-1+(j*numberOfNodes)] + 1;
                        }
                    //if width == linkWeights then...
                    }else{
                        //if (width or linkWeights (this is the same value -> watch condition above))>another width then...
                        if ((linkWeights[k+j*numberOfNodes] > width[i+k*numberOfNodes]) && (links[k+j*numberOfNodes] == 1)){
                            width[i+k*numberOfNodes] = linkWeights[k+j*numberOfNodes];
                            length[i+k*numberOfNodes] = length[i-1+(j*numberOfNodes)] + 1;
                        }
                    }

                }
            }
            pMax[i] = (pow(2, (gamma*length[(destination-1)*numberOfNodes+(i)])) - 1)/width[(destination-1)*numberOfNodes+(i)] ;
            //find min value for pMax and save it into array and save also her index
            if (pMax[i] < pMaxMin[exper]){
                pMaxMin[exper] = pMax[i];
                index[exper] = i;
            }
        }

        //toc
        toc = clock();
        measuredTime[exper] = toc - tic;
        //printf("\nmeasured time[%d] = %.6f", exper,  (double) measuredTime[exper]/CLOCKS_PER_SEC);
        directLinkPower[exper] = (pow(2, gamma)- 1)/linkWeights[(source-1)*numberOfNodes + (destination-1)];
    }

    printf("Average max power =              %.9f", sum(pMaxMin, totalExperiments)/totalExperiments);
    printf("\nAverage max power [dB] =         %.9f", 10*log10(sum(pMaxMin, totalExperiments)/totalExperiments));
    printf("\nAverage direct link power =      %.9f", sum(directLinkPower, totalExperiments)/totalExperiments);
    printf("\nAverage direct link power [dB] = %.9f", 10*log10(sum(directLinkPower, totalExperiments)/totalExperiments));
    //Average power savings = 100*(Average direct link power - Average max power)/Average direct link power
    printf("\nAverage power savings =          %.9f", 100*
        ((sum(directLinkPower, totalExperiments)/totalExperiments) - (sum(pMaxMin, totalExperiments)/totalExperiments))/
        (sum(directLinkPower, totalExperiments)/totalExperiments));
    printf("\nAverage time [s] =               %.3f", ((double) (timeSum(measuredTime, totalExperiments)/CLOCKS_PER_SEC))
        /totalExperiments);
    printf("\nAverage hop count =              %.3f", (double) sumInt(index, totalExperiments)/totalExperiments);

    //here should be two graphs:
    //plot(10*log10(P_max_min),'-') and
    //plot(10*log10(direct_link_power),':')


    return 0;
}
//return summation of all elements of integer array
int sumInt(int tab[], int numberOfElements){
    int i;
    double sum = 0;
    for(i = 0; i < numberOfElements; i++){
        sum = sum + tab[i];
    }
    return sum;
}

//return summation of all elements of time array
clock_t timeSum(clock_t tab[], int numberOfElements){
    int i;
    double sum = 0;
    for(i = 0; i < numberOfElements; i++){
        sum = sum + tab[i];
    }
    return sum;
}

//return summation of all elements of double array
double sum(double tab[], int numberOfElements){
    int i;
    double sum = 0;
    for(i = 0; i < numberOfElements; i++){
        sum = sum + tab[i];
    }
    return sum;
}

//compare and return smaller variable
double min(double x, double y){
    if (x < y){
        return x;
    }else{
        return y;
    }
}

//compare and return bigger variable
double max(double x, double y){
    if (x > y){
        return x;
    }else{
        return y;
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
