#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

typedef struct _particle
{
  double posX, posY, mass, velX, velY, b, Fx, Fy;
} particle;

typedef struct _particleBox
{
  struct _particleBox *nw;
  struct _particleBox *ne;
  struct _particleBox *sw;
  struct _particleBox *se;

  particle *star;

  double x;
  double y;
  double side;

  double mass;
  double centerOfMassX;
  double centerOfMassY;

  int isEmpty;
  int isLeaf;
} particleBox;

particle ** read_particle(int N, char * filename) {
  FILE* file = fopen(filename, "rb");
  double * buffer = (double*)malloc(sizeof(double)*6*N);

  int readOk = fread(buffer, sizeof(double)*6*N, 1, file);

  fclose(file);

  if(readOk != 1){
    printf("Error reading file");
    return NULL;
  }
  particle **array = (particle**)malloc(N * sizeof(particle*));


  for(int i = 0; i<N; i++){
    array[i] = (particle*)malloc(sizeof(particle));
  }
  int offset = 0;
  for(int i = 0; i<N; i++){
    array[i]->posX = *(buffer+offset);
    offset++;
    array[i]->posY = *(buffer+offset);
    offset++;
    array[i]->mass = *(buffer+offset);
    offset++;
    array[i]->velX = *(buffer+offset);
    offset++;
    array[i]->velY = *(buffer+offset);
    offset++;
    array[i]->b = *(buffer+offset);
    offset++;
    array[i]->Fx = 0;
    array[i]->Fy = 0;
  }

  free(buffer);
  buffer = NULL;
  return array;

}

void writeToFile(particle ** array, int N) {
  char * filename = "result.gal";
  FILE * fp;
  fp = fopen(filename, "wb");
  for (int i = 0; i < N; i++) {
    fwrite(array[i], sizeof(double)*6, 1, fp);
  }
  fclose(fp);
}


particleBox * createBox(double x, double y, double side){
  particleBox *box;

  box = malloc(sizeof(particleBox));

  box->mass = 0;
  box->centerOfMassX = 0;
  box->centerOfMassY = 0;

  box->x = x;
  box->y = y;
  box->side = side;

  box->star = NULL;
  box->isEmpty = 1;
  box->isLeaf = 1;


  return box;
}

void fitParticle(particleBox * box, particle * star){
  if(box->isEmpty == 1){
    //box->mass += star->mass;
    //printf("Mass for box %f, %f: %f\n",box->x, box->y, box->mass);

    box->star = star;
    box->isEmpty = 0;
    box->isLeaf = 1;
  }
  else if(box->isLeaf == 1){
    particle *oldStar = box->star;

    box->nw = createBox(box->x-box->side/4, box->y-box->side/4, box->side/2);
    box->ne = createBox(box->x+box->side/4, box->y-box->side/4, box->side/2);
    box->sw = createBox(box->x-box->side/4, box->y+box->side/4, box->side/2);    
    box->se = createBox(box->x+box->side/4, box->y+box->side/4, box->side/2);
    box->star = NULL;
    box->isLeaf = 0;
    fitParticle(box, oldStar);
    fitParticle(box, star);      
  }else{
    if(star->posX < box->x && star->posY < box->y){
      fitParticle(box->nw, star);
      }else if(star->posX >= box->x && star->posY < box->y){
        fitParticle(box->ne, star);
      }else if(star->posX < box->x && star->posY >= box->y){
        fitParticle(box->sw, star);
      }else if(star->posX >= box->x && star->posY >= box->y){
        fitParticle(box->se, star);
      }
    else {
      printf("Error when trying to insert particle");
    }

  }
}


void calcCenterOfMass(particleBox * box){
  double r1, r2, r3, r4;
  double m1, m2, m3, m4;
  double xCenter, yCenter;

  r1 = box->nw->centerOfMassX*box->nw->mass;
  r2 = box->ne->centerOfMassX*box->ne->mass;
  r3 = box->sw->centerOfMassX*box->sw->mass;
  r4 = box->se->centerOfMassX*box->se->mass;

  m1 = box->nw->mass;
  m2 = box->ne->mass;
  m3 = box->sw->mass;
  m4 = box->se->mass; 

  xCenter = (r1+ r2 + r3 + r4)/(m1 + m2 + m3 + m4);

  r1 = box->nw->centerOfMassY;
  r2 = box->ne->centerOfMassY;
  r3 = box->sw->centerOfMassY;
  r4 = box->se->centerOfMassY;

  yCenter = (r1*m1 + r2*m2 + r3*m3 + r4*m4)/(m1 + m2 + m3 + m4);

  box->centerOfMassX = xCenter;
  box->centerOfMassY = yCenter;
}

void calcMass(particleBox * box){
  if(box->isEmpty){
    /* box->mass = 0;
    box->centerOfMassX = 0;
    box->centerOfMassY = 0; */

  }else if(box->isLeaf){
    box->mass=box->star->mass;
    box->centerOfMassX = box->star->posX;
    box->centerOfMassY = box->star->posY;

  }else{
    calcMass(box->nw);
    calcMass(box->ne);
    calcMass(box->sw);
    calcMass(box->se);

    double mass = box->nw->mass + box->ne->mass + box->sw->mass + box->se->mass; 
    box->mass += mass;

    calcCenterOfMass(box);
  }
}

void calcForce(particle *star, particleBox *box, double thetaMax){
  double rx, ry, denom;
  double e0 = 0.001;
  rx = star->posX - box->centerOfMassX;
  ry = star->posY - box->centerOfMassY;
  double r = sqrt(rx*rx+ry*ry);
 // if(box->star != star){
/*     if(box->isEmpty){
      
    } */
     if(box->isLeaf || (box->side / r) <= thetaMax){
      denom = (r+e0)*(r+e0)*(r+e0);
      star->Fx += box->mass * rx / denom;
      star->Fy += box->mass * ry / denom;
    } 
    else {
      calcForce(star, box->nw, thetaMax);
      calcForce(star, box->ne, thetaMax);
      calcForce(star, box->sw, thetaMax);
      calcForce(star, box->se, thetaMax);
    }
  //}
}

void deleteBoxes(particleBox * box) {
  if (box->isLeaf) {
    free(box);
    box = NULL;
  } else {
    deleteBoxes(box->nw);
    deleteBoxes(box->ne);
    deleteBoxes(box->sw);
    deleteBoxes(box->se);
    free(box);
    box = NULL;
  }
}

void printArray(particle ** a, int N) {
  for (int i = 0; i < N; i++) {
    printf("mass: %f, Xpos: %f, Ypos: %f, b: %f \n", a[i]->mass, a[i]->posX, a[i]->posY, a[i]->b);
  }
}

particleBox* buildTree(particle ** array, int N) {
    particleBox * root;
    root = malloc(sizeof(particleBox));
  
    root->mass = 0;
    root->x = 0.5;
    root->y = 0.5;
    root->side = 1;
    root->star = NULL;

    root->nw = NULL;
    root->ne = NULL;
    root->sw = NULL;
    root->se = NULL;

    root->isEmpty = 1;
    root->isLeaf = 1;

    root->centerOfMassX = 0.5;
    root->centerOfMassY = 0.5;
    for(int i = 0; i<N; i++){
      fitParticle(root, array[i]);
    }
    return root;
}

int main(int argc, char* argv[]){  
  if (argc != 8){
      printf("Wrong number of input arguments\n");
      return 1;
  }  
  int N = atoi(argv[1]);
  char* filename = argv[2];
  int n_steps = atoi(argv[3]);
  double delta_t = atof(argv[4]);
  double theta_max = atof(argv[5]);
  int graphics = atoi(argv[6]);
  int n_threads = atoi(argv[7]);
  printf("Command line arguments given: %d, %s, %d, %f, %f, %d, %d \n", N, filename, n_steps, delta_t, theta_max, graphics, n_threads);
  const double G = 100.0 / N;

  omp_set_num_threads(n_threads);

  particle **array = read_particle(N, filename);
  //printArray(array, N);
  particleBox *root = NULL;
  for(int i = 0; i<n_steps; i++) {
    root = buildTree(array, N);
    calcMass(root);

    #pragma omp parallel for schedule(static)
    for(int j =0; j<N; j++){
      array[j]->Fx = 0;
      array[j]->Fy = 0;

      calcForce(array[j], root, theta_max);

      array[j]->Fx *= -G;
      array[j]->Fy *= -G;
    }
    deleteBoxes(root);
    //#pragma omp parallel for schedule(static) num_threads(n_threads)
    for (int j = 0; j < N; j++) {
      array[j]->velX += delta_t*array[j]->Fx;
      array[j]->velY += delta_t*array[j]->Fy;
      array[j]->posX += delta_t*array[j]->velX;
      array[j]->posY += delta_t*array[j]->velY;
    }
  }
  //printf("\n");
  //printArray(array, N);
  writeToFile(array, N);
  for (int i = 0; i < N; i++) {
    free(array[i]);
    array[i] = NULL;
  }
  free(array);
  array = NULL;
  return 0;
}

//./galsim 4 input_data/circles_N_4.gal 500 1e-5 0.1 0
//./galsim 10 input_data/ellipse_N_00010.gal 200 1e-5 0.1 0

//./galsim 10 input_data/ellipse_N_00010.gal 5 1e-5 0.1 0

//./compare_gal_files 10 ../result.gal ../ref_output_data/ellipse_N_00010_after200steps.gal