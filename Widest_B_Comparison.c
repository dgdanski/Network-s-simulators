//author: Daniel Gdanski
//Sharjah, UAE, 15.11.2016

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

double min(double x, double y);
int sum(int labeled[], int numberOfNodes);
double randn (double mu, double sigma);

int main(){
    int numberOfNodes = 20;
    //source nodes
    int source = 1;
    //destination nodes
    int destination = 4;
    //infinity
    int inf = 2147483647;
    double x[numberOfNodes], y[numberOfNodes];
    int powerOfNumberOfNodes = numberOfNodes * numberOfNodes, orderNumber = 1,
        numberOfZeros[powerOfNumberOfNodes/(numberOfNodes+1)], cost[powerOfNumberOfNodes];
    int i, j, k, m, z, tt, counter, index, labeled[numberOfNodes], rPath[powerOfNumberOfNodes];
    double c, maxShadow, rWidest, maxR, temp, minimumLinkWeight, maxRl, referenceDistance = 100,
        width[powerOfNumberOfNodes], length[powerOfNumberOfNodes], linkWeights[powerOfNumberOfNodes],
        linkCosts[powerOfNumberOfNodes], alpha[powerOfNumberOfNodes], originalAlpha[powerOfNumberOfNodes],
        euclid[numberOfNodes * numberOfNodes], shadow[powerOfNumberOfNodes], r[numberOfNodes], dist[numberOfNodes],
        prev[numberOfNodes], rl[powerOfNumberOfNodes];
    bool links[powerOfNumberOfNodes];
    clock_t beginTimeWidest, endTimeWidest, beginTimeB, endTimeB;

    if (numberOfNodes < destination){
        printf("Number of nodes can't be less than destination!\nChange the value and please try again.\n\n");
        return 0;
    }

    srand(time(0));
    //Generating the nodes, assume 20 nodes, randomly scattered in a 100*100 area
    for(i = 0; i < numberOfNodes; i++){
        x[i] = rand()%101 + ((double)rand()/(double)RAND_MAX);
        y[i] = rand()%101 + ((double)rand()/(double)RAND_MAX);
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

    c = 1/maxShadow;

    //generating the random link weights
    for(i = 0; i < powerOfNumberOfNodes; i++){
        if (links[i] == 1){
            linkWeights[i] = log2(1+(pow(10,4))*c*shadow[i]*pow((euclid[i]/referenceDistance), -4));
            alpha[i] = linkWeights[i];
            originalAlpha[i] = linkWeights[i];
        }else{
            linkWeights[i] = 0;
            originalAlpha[i] = 0;
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

    //fill rl with 0
    for(i=0; i<powerOfNumberOfNodes; i++){
        rl[i] = 0;
    }

    //tic
    beginTimeWidest = clock();

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
        r[i] = width[(destination-1)*numberOfNodes+(i)]/length[(destination-1)*numberOfNodes+(i)];
    }

    //toc
    endTimeWidest = clock();
    printf("Time widest = %.6f seconds", (double)((endTimeWidest - beginTimeWidest)/ CLOCKS_PER_SEC));

    //find max(r)
    maxR = r[0];
    for(i = 0; i < numberOfNodes; i++){
        if (r[i] > maxR){
            maxR = r[i];
            index = i;
        }
    }
    rWidest = maxR;
    printf("\nR_widest = %.6f", rWidest);

    //Comparing against Algorithm B
    //tic
    beginTimeB = clock();
    counter = 0;
    for(i = 0; i < powerOfNumberOfNodes; i++){
       if (originalAlpha[i] > originalAlpha[(source-1)*numberOfNodes+(destination-1)]){
            counter++;

            //clear linkCosts and linkWeights
            for(j = 0; j < powerOfNumberOfNodes; j++){
                linkCosts[j] = 0;
                linkWeights[j] = 0;
            }


            j = 0;
            orderNumber = 1;
            for(k = 0; k < powerOfNumberOfNodes; k++){
                if(orderNumber == numberOfZeros[j]){
                    links[k] = 0;
                    cost[k] = 0;
                    j++;
                }else{
                    links[k] = 1;
                    cost[k] = 1;
                }
                orderNumber++;
                alpha[k] = originalAlpha[k];
            }

            //Remove all links for which log(1+PG_{l}/N_{0}B) < a
            for(j = 0; j < powerOfNumberOfNodes; j++){
                if(alpha[j] < alpha[i]){
                    links[j] = 0;
                    cost[j] = inf;
                }
            }

            //labeled[j] = 1 indicates that node i is permanently labeled, initially
            //labeled[source - 1] = 1 only (source - 1 is the only permanently labeled node)
            for(j = 0; j < numberOfNodes; j++){
                labeled[j] = 0;
            }
            labeled[source - 1] = 1;

            //Dist[j] is the distance from the source to node j (found so far). The following are the initial distances
            for(j = 0; j < numberOfNodes; j++){
                dist[j] = (double) cost[(source-1)*numberOfNodes + j];
            }



            //Prev[j] indicates the node previous to j along the shortest path from source (to j)
            //Initially, Prev[j] = source for all nodes j that are neighbors of source
            for(j = 0; j < numberOfNodes; j++){
                if(cost[(source-1)*numberOfNodes + j] < inf && cost[(source-1)*numberOfNodes + j] > 0){
                    prev[j] = source;
                }else{
                    prev[j] = 0;
                }
            }
            z = 0;
            while(sum(labeled, numberOfNodes) < numberOfNodes){
                index = -1;
                //Finding the next permanently labeled node
                temp = inf;
                for(j = 0; j < numberOfNodes; j++){
                    if(labeled[j] == 0){
                        if (dist[j] < temp){
                            temp = dist[j];
                            index = j;
                        }
                    }
                }

                if(index != -1){
                    labeled[index] = 1;
                }else{
                    break;
                }

                //Updating labels
                for(j = 0; j < numberOfNodes; j++){
                    //updating total cost of shortest path
                    if(labeled[j] == 0){
                        if (dist[index] + cost[(index * numberOfNodes) + j] < dist[j]){
                            prev[j] = index + 1;
                        }
                        dist[j] = min(dist[j], dist[index] + (double) cost[(index * numberOfNodes) + j]);
                    }
                }
                z++;
            }

            //Construct actual path (sequence of links)
            if (dist[destination - 1] < inf){
                rPath[0] = destination;
                tt = destination;
                j = 1;
                while (source != tt){
                    rPath[j] = prev[tt - 1];
                    tt = prev[tt - 1];
                    j++;
                }

                //Recover costs and weights of the links that constitute to the shortest path
                minimumLinkWeight = inf;
                for(k = 0; k < j - 1; k++){
                    linkCosts[k] = cost[(rPath[k+1]-1)*numberOfNodes + rPath[k]-1];
                    linkWeights[k] = alpha[(rPath[k+1]-1)*numberOfNodes + rPath[k]-1];
                    //Display costs and weights of the individual links of the shortest path
                    //printf("\nCost   of link [%d] = %.2f", k+1, linkCosts[k]);
                    //printf("\nWeight of link [%d] = %.2f", k+1, linkWeights[k]);
                    if (linkWeights[k] < minimumLinkWeight){
                        minimumLinkWeight = linkWeights[k];
                    }
                }
                rl[counter] = minimumLinkWeight/dist[destination-1];
            }else{
                rl[counter] = 0;
            }
        }
    }

    //toc
    endTimeB = clock();
    printf("\nTime B = %.6f", (double)(((endTimeB - beginTimeB)/ CLOCKS_PER_SEC)));

    counter++;
    rl[counter] = originalAlpha[(source-1)*numberOfNodes + (destination - 1)];
    maxRl = 0;
    for(i = 1; i < counter+1; i++){
        if(rl[i] > maxRl){
            maxRl = rl[i];
        }
    }
    printf("\nR_B = %.6f", maxRl);

    if ((endTimeB - beginTimeB) != 0){
        printf("\nTime ratio = %.9f", (double)(((endTimeB - beginTimeB) - (endTimeWidest - beginTimeWidest))/(endTimeB - beginTimeB)));
    }
    printf("\nSpectral efficiency ratio = %.2f", rWidest/maxRl);

    return 0;
}

//compare and return smaller variable
double min(double x, double y){
    if (x < y){
        return x;
    }else{
        return y;
    }
}

//return summation of all elements of an array
int sum(int labeled[], int numberOfNodes){
    int i;
    double sum = 0;
    for(i = 0; i < numberOfNodes; i++){
        sum = sum + labeled[i];
    }
    return sum;
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
