#include <BasicLinearAlgebra.h>
#include <Servo.h>
#include <math.h>
using namespace BLA;
Servo joint1;
Servo joint2;

int joint1Pin = 5;
int joint2Pin = 6;

Matrix<3,3> Jcbn;

Matrix<3,3> JcbnT;
Matrix<3,3> Iden;
Matrix<3,3> damp;
Matrix<3,3> Invdamp;

Matrix<3,1> DeltaP;
Matrix<3,1> Delta0;


double a1;
double a2; //measured from the horizontal axis touching the servo
double b;
double xd; //x distance
double yd;
double zd;

int L1 = 80;//mm
int L2 = 80;//mm
int P = 1;

//home position
float xii = 54.72;
float yii = 0.00;
float zii = 0.00;

int stepCount; //stepper motor step counter

unsigned long lastTime = 0;
double dt = 0;
Matrix<3,1> iP = {0,0,0};

Matrix<3,1> TargetP = {0,0,0};
Matrix<3,1> Position = {0,0,0};
Matrix<3,1> DirectionS = {0,0,0};

float sigmaMin = 1;
float s = 40;
float t = 0;
Matrix<3,1> ErrorP = {0,0,0};
Matrix<3,1> getDirection(Matrix<3,1> DirectionS, Matrix<3,1> Position, Matrix<3,1>* ErrorP, Matrix<3,1>* TargetP, Matrix<3,1>* iP, double dt, float* sigmaMin);

void setup()
{
  Serial.begin(9600);
  joint1.attach(joint1Pin);
  joint1.write(70);
  joint2.attach(joint2Pin);
  joint2.write(40);
  pinMode(8,OUTPUT); //stepper motor pins
  pinMode(9,OUTPUT);
  pinMode(10,OUTPUT);
  pinMode(11,OUTPUT);
  pinMode(joint1Pin, OUTPUT);//a1
  pinMode(joint2Pin, OUTPUT);//a2
  pinMode(51,INPUT); //button pin
  pinMode(50,INPUT); //also button pin
  pinMode(49,INPUT); //button pin
  pinMode(48,INPUT); //also button pin
  pinMode(47,INPUT); //button pin
  pinMode(46,INPUT); //also button pin
  pinMode(45,INPUT);

  stepCount = 0;
  b = 0;
  a1 = joint1.read();
  a2 = joint2.read()-180+a1;
  
  FK(L1,L2,a1,a2,b,&xd,&yd,&zd);
  Position = {xd,yd,zd};
  Serial.println(Position);
  
}

void loop()
{
  unsigned long now = micros();
  dt = (now - lastTime) * 1e-6;
  lastTime = now;

  FK(L1,L2,a1,a2,b,&xd,&yd,&zd);
  Position = {xd,yd,zd};
  Matrix<3,1> Angle = {a1,a2,b};
  
  if(digitalRead(45) == HIGH){
    FK(L1,L2,a1,a2,b,&xd,&yd,&zd);
    Position = {xd,yd,zd};
    Serial.print("Position:"); Serial.println(Position);
    Matrix<3,1> xError = {0,yd-yii,zd-zii};
    Serial.print("xError:"); Serial.println(xError);
    Matrix<3,1> yError = {xd-xii,0,zd-zii};
    Serial.print("yError:"); Serial.println(yError);
    Matrix<3,1> zError = {xd-xii,yd-yii,0};
    Serial.print("zError:"); Serial.println(zError);


  }
  if(digitalRead(51) == HIGH)
  {
    DirectionS = {1*s, //x
                  0*s, //y
                  0*s};//z

    b = IKnew(L1,L2,a1,a2,b,&stepCount,xd,yd,zd,DirectionS,&sigmaMin);
    a1 = joint1.read();
    a2 = joint2.read()-180+a1;
    
    delay(1);
  }

  if(digitalRead(50) == HIGH)
  {
    DirectionS = {-1*s, //x
                  0*s, //y
                  0*s};//z
    b = IKnew(L1,L2,a1,a2,b,&stepCount,xd,yd,zd,DirectionS,&sigmaMin);
    a1 = joint1.read();
    a2 = joint2.read()-180+a1;
    
    delay(1);
  }

  if(digitalRead(49) == HIGH)
  {
    
    DirectionS = {0*s, //x
                  1*s, //y
                  0*s};//z
  
    b = IKnew(L1,L2,a1,a2,b,&stepCount,xd,yd,zd,DirectionS,&sigmaMin);
    a1 = joint1.read();
    a2 = joint2.read()-180+a1;
    
    
    delay(1);
  }

  if(digitalRead(48) == HIGH)
  {
    DirectionS = {0*s, //x
                  -1*s, //y
                  0*s};//z
    b = IKnew(L1,L2,a1,a2,b,&stepCount,xd,yd,zd,DirectionS,&sigmaMin);
    a1 = joint1.read();
    a2 = joint2.read()-180+a1;
    
    delay(1);
  }

  if(digitalRead(47) == HIGH)
  {
    DirectionS = {0*s, //x
                  0*s, //y
                  1*s};//z
    b = IKnew(L1,L2,a1,a2,b,&stepCount,xd,yd,zd,DirectionS,&sigmaMin);
    a1 = joint1.read();
    a2 = joint2.read()-180+a1;
    
    delay(1);
  }

  if(digitalRead(46) == HIGH)
  { 
    DirectionS = {0*s, //x
                  0*s, //y
                  -1*s};//z
    b = IKnew(L1,L2,a1,a2,b,&stepCount,xd,yd,zd,DirectionS,&sigmaMin);
    a1 = joint1.read();
    a2 = joint2.read()-180+a1;
    
    delay(1);
  }

  
  
  
}

double IKnew(int L1, int L2, double a1, double a2, double b, int* stepCount, double xi, double yi, double zi, Matrix<3,1> DirectionF, float* sigmaMin){
  double lambda,xd,yd,zd,xf,yf,zf;
  if(DirectionF(0,0) == 0 && DirectionF(1,0) == 0 && DirectionF(2,0) == 0){ //check for no movement asked
    return b;
  }

  b = b*(PI/180);
  a1 = a1*(PI/180);
  a2 = a2*(PI/180);

  //Jacobian row x column;
  double J11 = -L1*sin(a1)*cos(b);
  double J12 = -L2*sin(a2)*cos(b);
  double J13 = -sin(b)*(L1*cos(a1)+L2*cos(a2));
  double J21 = -L1*sin(a1)*sin(b);
  double J22 = -L2*sin(a2)*sin(b);
  double J23 = cos(b)*(L1*cos(a1)+L2*cos(a2));
  double J31 = L1*cos(a1);
  double J32 = L2*cos(a2);
  double J33 = 0;

  Jcbn = {J11, J12, J13,
          J21, J22, J23,
          J31, J32, J33};
  JcbnT = {J11, J21, J31,
           J12, J22, J32,
           J13, J23, J33};

  
  *sigmaMin = checkSingularity(Jcbn, JcbnT);
  if(*sigmaMin < 0.1){ //emergency home
    joint1.write(70);
    joint2.write(40);
    Serial.print("return to home");
    return 0;
  }  

  lambda = 20 + (2000/(((*sigmaMin)*(*sigmaMin))+ 1e-1)); //adjust damping coefficient
  
  Iden = {lambda*lambda, 0, 0,
          0, lambda*lambda, 0,
          0, 0, lambda*lambda};
  

  Matrix<3,3> JtJ = JcbnT*Jcbn;
  damp = JtJ + Iden; //Jacobian*Jacobian Transpose + Lambda^2 * Identity Matrix

  Invdamp = Inverse(damp); 
  



  
  

  DeltaP = getDirection(DirectionF, Position, &ErrorP, &TargetP, &iP, dt, sigmaMin); //passes directionF into PID controller
  DeltaP = DeltaP + (DeltaP * Matrix<1,1> {-1/(*sigmaMin+0.1)}); //adjusts change in position according to proximity to singular regions
  if((abs(DeltaP(2,0) + zi)) > (abs(150*sin(a1)))){ //detects the outer boundary
    if((abs(DeltaP(0,0) + xi) <= abs(1.2*xi)) && (abs(DeltaP(1,0) + yi) <= abs(1.2*yi)) || (abs(DeltaP(2,0) + zi)) <= (abs(1.2*zi))){//detects if user input will move end effector into our away from outer singular region

    }
    else{
    Serial.print("Too close to edge");
    DeltaP = {0,0,0};
    }

  }
  else if((abs(DeltaP(0,0) + xi)) < 8 && (abs(DeltaP(1,0) + yi)) < 8){//detects if user input will move end effector into our away from central cylindrical singularity
    if((DeltaP(0,0) + xi) >= xi && (DeltaP(1,0) + yi >= yi)){

    }
    else{
    Serial.print("Too close to center");
    DeltaP = {0,0,0};
    }
  }


  Delta0 = Invdamp*JcbnT*DeltaP;  //compute change in angles required
  Matrix<3,1> ActualDeltaP = Jcbn * Delta0; //compute what the change in position shouldve been
  Matrix<3,1> Residual = DeltaP - ActualDeltaP; //gather residual error from DLS
  Matrix<3,1> Delta0Correction = Invdamp * JcbnT * Residual; 
  Delta0 = Delta0 + Delta0Correction; 
  //resulting error can be corrected by PID controller without fighting against the jacobian as its error has already been accounted for
  
  float a1F = Delta0(0,0) + a1;
  float a2F = Delta0(1,0) + a2;
  float bF = Delta0(2,0) + b;


  bF = bF*(180/PI); //convert back into degrees
  a1F = a1F*(180/PI);
  a2F = a2F*(180/PI);




  if(a1F-int(a1F) > 0.5){ //rounding formula since my servos only turn to an integer angle
    joint1.write(a1F+1);
  }
  else{
    joint1.write(a1F);
  }
  if(a2F-int(a2F) > 0.5){
    joint2.write(180-(a1F-a2F)+1);
  }
  else{
    joint2.write(180-(a1F-a2F));
  }
  if(abs(Delta0(2,0)) <= 0.0002){
    return(bF);
  }
  if(Delta0(2,0) > 0){ //this is true if is positive
    if(bF-int(bF) > 0.5){
      //clockwise(bF);
      *stepCount = *stepCount + clockdumb(Delta0(2,0), &P)+1; //calls counterclockwise rotation with respect to the precalculated change in b angle
    }
    else{
      *stepCount = *stepCount + clockdumb(Delta0(2,0), &P);
    }
  }
  else{ //this is true if bF is negative
    if(abs(bF-int(bF)) > 0.5){ 
      
      *stepCount = *stepCount - clockwise(abs(Delta0(2,0)), &P)+1; //same but clockwise
    }
    else{
      *stepCount = *stepCount - clockwise(abs(Delta0(2,0)), &P);
    }
  }
  //Serial.println(*stepCount);
  return(bF);
}

int clockdumb(float change, int* P){ //counterclockwise spin to the change in b angle

//4096 steps per rotation
    int e = 1;
    int step = *P;
    double runs = change*(180/PI)*(4096/360); 
    int localStepCount = 0;
    if(runs - int(runs) > 0.5){ //rounds into integer
      runs = int(runs)+1;
    }
    else{
      runs = int(runs);
    }
  while(localStepCount <= runs){
    switch(step){ //each case activates in order and is counted and added to the local step count
      case 1:
      {
        digitalWrite(8, HIGH);
        digitalWrite(9, LOW);
        digitalWrite(10, LOW);
        digitalWrite(11, HIGH);
        delay(e);
        step++;
        localStepCount++;
        *P = 1;
        
        break;
      }
      case 2:
      {
        digitalWrite(8, LOW);
        digitalWrite(9, LOW);
        digitalWrite(10, LOW);
        digitalWrite(11, HIGH);
        delay(e*2);
        step++;
        localStepCount++;
        *P = 2;
        break;
      }
      case 3:
      {
        digitalWrite(8, LOW);
        digitalWrite(9, LOW);
        digitalWrite(10, HIGH);
        digitalWrite(11, HIGH);
        delay(e);
        step++;
        localStepCount++;
        *P=3;
        break;
      }
      case 4:
      {
        digitalWrite(8, LOW);
        digitalWrite(9, LOW);
        digitalWrite(10, HIGH);
        digitalWrite(11, LOW);
        delay(e*2);
        step++;
        localStepCount++;
        *P=4;
        break;
      }
      case 5:
      {
        digitalWrite(8, LOW);
        digitalWrite(9, HIGH);
        digitalWrite(10, HIGH);
        digitalWrite(11, LOW);
        delay(e);
        step++;
        localStepCount++;
        *P=5;
        break;
      }
      case 6:
      {
        digitalWrite(8, LOW);
        digitalWrite(9, HIGH);
        digitalWrite(10, LOW);
        digitalWrite(11, LOW);
        delay(e*2);
        step++;
        localStepCount++;
        *P=6;
        break;
      }
      case 7:
      {
        digitalWrite(8, HIGH);
        digitalWrite(9, HIGH);
        digitalWrite(10, LOW);
        digitalWrite(11, LOW);
        delay(e);
        step++;
        localStepCount++;
        *P=7;
        break;
      }
      case 8:
      {
        digitalWrite(8, HIGH);
        digitalWrite(9, LOW);
        digitalWrite(10, LOW);
        digitalWrite(11, LOW);
        delay(e*2);
        step++;
        localStepCount++;
        *P=8;
        break;
      }
      default: //resets the pattern should the change in angle demand it
      {
        step = step%8;
        break;

      }

    }
  }
    delay(e);
    return localStepCount; //returns how many steps were traveresed in this instance of the function
    
}




int clockwise(float change, int* P){ //same as previous for a clockwise rotation
    int e = 1;
    int step = *P;
    
    float runs = change*(180/PI)*(4096/360); 
    
    int localStepCount = 0;
    if(runs - int(runs) > 0.5){
      runs = int(runs)+1;
    }
    else{
      runs = int(runs);
    }
    
  while(localStepCount <= runs){ 
    switch(step){
      case 1:
      {
        digitalWrite(8, HIGH);
        digitalWrite(9, LOW);
        digitalWrite(10, LOW);
        digitalWrite(11, LOW);
        delay(e*2);
        step++;
        localStepCount++;
        *P=1;
        break;
      }
      case 2:
      {
        digitalWrite(8, HIGH);
        digitalWrite(9, HIGH);
        digitalWrite(10, LOW);
        digitalWrite(11, LOW);
        delay(e);
        step++;
        localStepCount++;
        *P=2;
        break;
      }
      case 3:
      {
        digitalWrite(8, LOW);
        digitalWrite(9, HIGH);
        digitalWrite(10, LOW);
        digitalWrite(11, LOW);
        delay(e*2);
        step++;
        localStepCount++;
        *P=3;
        break;
      }
      case 4:
      {
        digitalWrite(8, LOW);
        digitalWrite(9, HIGH);
        digitalWrite(10, HIGH);
        digitalWrite(11, LOW);
        delay(e);
        step++;
        localStepCount++;
        *P=4;
        break;
      }
      case 5:
      {
        digitalWrite(8, LOW);
        digitalWrite(9, LOW);
        digitalWrite(10, HIGH);
        digitalWrite(11, LOW);
        delay(e*2);
        step++;
        localStepCount++;
        *P=5;
        break;
      }
      case 6:
      {
        digitalWrite(8, LOW);
        digitalWrite(9, LOW);
        digitalWrite(10, HIGH);
        digitalWrite(11, HIGH);
        delay(e);
        step++;
        localStepCount++;
        *P=6;
        break;
      }
      case 7:
      {
        digitalWrite(8, LOW);
        digitalWrite(9, LOW);
        digitalWrite(10, LOW);
        digitalWrite(11, HIGH);
        delay(e*2);
        step++;
        localStepCount++;
        *P=7;
        break;
      }
      case 8:
      {
        digitalWrite(8, HIGH);
        digitalWrite(9, LOW);
        digitalWrite(10, LOW);
        digitalWrite(11, HIGH);
        delay(e);
        step++;
        localStepCount++;
        *P=8;
        break;
      }
      default:
      {
        step = step%8;
        break;

      }

    }
  }
    delay(e);
    return localStepCount;
    
}
void FK(int L1, int L2, double a1, double a2, double b, double* xd, double* yd, double* zd)
{
  double c1,c2,c3,s1,s2,s3;
  b = b*(PI/180);
  a1 = a1*(PI/180);
  a2 = a2*(PI/180);
  
  //preloading trigonometric functions
  c1 = cos(b); s1 = sin(b);
  c2 = cos(a1); s2 = sin(a1);
  c3 = cos(a2); s3 = sin(a2);

      
  *xd = L2*c3*c1 + L1*c2*c1; 
  *yd = L2*c3*s1 + L1*c2*s1;
  *zd = L2*s3     + L1*s2;

}

float checkSingularity(Matrix <3,3> Jcbn, Matrix <3,3> JcbnT)
{
  //the purpose of this function is to find the smallest singular value of the JJT, in SVD the smallest singular value correlates to the 
  //smallest U and V hence the smallest shift in position
  Matrix<3,3> A = Jcbn*JcbnT;
  Matrix<3,1> w;
 //solve Ax = w where x is any random vector
  
  Matrix<3,1> x = {1,1,1};
  Matrix<3,3> IdenSmall = {0.00001, 0, 0,
                          0, 0.00001, 0,
                          0, 0, 0.00001};
  // Regularization to avoid singularity
  A = A + IdenSmall;
  Matrix<3,3> Ai = Inverse(A);
  for(int i = 0; i<4; i++){
    // Solve A w = x using the inverse
  
    w = Ai * x;
    float mag = 1/sqrt(w(0,0)*w(0,0) + w(1,0)*w(1,0) + w(2,0)*w(2,0));
    w = {w(0,0)*mag, w(1,0)*mag, w(2,0)*mag};
    x = w;
  }
  
  Matrix<1,1> Eig = (~w)*(Inverse(A))*w;
  float sigmaMin = sqrt(Eig(0,0));

  Serial.print("sigmaMin:"); Serial.println(sigmaMin);//showing user proximity to singular region
  return sigmaMin;
  

}

Matrix<3,1> getDirection(Matrix<3,1> DirectionS, Matrix<3,1> Position, Matrix<3,1>* ErrorP, Matrix<3,1>* TargetP, Matrix<3,1>* iP, double dt, float* sigmaMin)
{ 
  Matrix<3,1> Target, TargetPrev;
  float kP = 0.06;
  float kD = 0.003;
  float kI = 0.0001;
  if(*sigmaMin < 10){ //near singular
    *ErrorP = 0;
    kD = 0;
    kI = 0;
    kP = 0.01;
  }
  if(*sigmaMin < 0.1){ //inside singular
    *ErrorP = 0;
    *TargetP = {0,0,0};
    kD = 0;
    kI = 0;
    kP = 0.0;
    *iP = 0;
    return {0,0,0};
  }
  


  TargetPrev = *TargetP; //TargetPrev used for simplification
  if(TargetPrev(0,0) == 0 && TargetPrev(1,0) == 0 && TargetPrev(2,0) == 0)
  {
    Target = Position;
  }
  else{

  Target = *TargetP + DirectionS*(Matrix<1,1> {dt}); //integrate directionS to get units of position

  }
  Matrix<3,1> actual = Position;
  
  Matrix<3,1> e = Target - actual;
  Matrix<3,1> d = (e - *ErrorP)*(Matrix<1,1> {1/dt});
  Matrix<3,1> i = *iP + e*(Matrix<1,1> {dt});
  
  Matrix<3,1> u = e*(Matrix<1,1> {kP}) + d*(Matrix<1,1> {kD}) + i*(Matrix<1,1> {kI});
  
  
  
  *ErrorP = e; //establish previous values
  *TargetP = Target;
  *iP = i;
  
  return u; //return correct velocity


}
