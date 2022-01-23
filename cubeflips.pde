// Processing code by Etienne JACOB
// motion blur template by beesandbombs
// opensimplexnoise code by Kurt Spencer is used
// --> code here : https://gist.github.com/Bleuje/fce86ef35b66c4a2b6a469b27163591e

int[][] result;
float t, c;

float c01(float x)
{
  return constrain(x,0,1);
}

float ease(float p) {
  return 3*p*p - 2*p*p*p;
}

float ease(float p, float g) {
  if (p < 0.5) 
    return 0.5 * pow(2*p, g);
  else
    return 1 - 0.5 * pow(2*(1 - p), g);
}

float map(float x, float a, float b, float c, float d, boolean constr)
{
  return constr ? c01(map(x,a,b,c,d)) : map(x,a,b,c,d);
}

float mp01(float x, float a, float b)
{
  return map(x,a,b,0,1,true);
}

float pow_(float p,float g)
{
  return 1-pow(1-p,g);
}

float tanh(float x)
{
  return (float)Math.tanh(x);
}

float softplus(float q,float p){
  float qq = q+p;
  if(qq<=0){
    return 0;
  }
  if(qq>=2*p){
    return qq-p;
  }
  return 1/(4*p)*qq*qq;
}

float mn = .5*sqrt(3), ia = atan(sqrt(.5));

void push() {
  pushMatrix();
  pushStyle();
}

void pop() {
  popStyle();
  popMatrix();
}

void draw() {

  if (!recording) {
    t = (mouseX*1.25/width)%1;
    c = mouseY*1.0/height;
    if (mousePressed)
      println(c);
    draw_();
  } else {
    for (int i=0; i<width*height; i++)
      for (int a=0; a<3; a++)
        result[i][a] = 0;

    c = 0;
    for (int sa=0; sa<samplesPerFrame; sa++) {
      t = map(frameCount-1 + sa*shutterAngle/samplesPerFrame, 0, numFrames, 0, 1);
      draw_();
      loadPixels();
      for (int i=0; i<pixels.length; i++) {
        result[i][0] += pixels[i] >> 16 & 0xff;
        result[i][1] += pixels[i] >> 8 & 0xff;
        result[i][2] += pixels[i] & 0xff;
      }
    }

    loadPixels();
    for (int i=0; i<pixels.length; i++)
      pixels[i] = 0xff << 24 | 
        int(result[i][0]*1.0/samplesPerFrame) << 16 | 
        int(result[i][1]*1.0/samplesPerFrame) << 8 | 
        int(result[i][2]*1.0/samplesPerFrame);
    updatePixels();
    
    if (frameCount<=numFrames)
    {
      saveFrame("fr###.gif");
      println(frameCount,"/",numFrames);
    }
    
    if (frameCount==numFrames)
      stop();
  }
}

//////////////////////////////////////////////////////////////////////////////

int samplesPerFrame = 7;
int numFrames = 440;        
float shutterAngle = .7;

boolean recording = false;

OpenSimplexNoise noise;

//int SEED = 1234;
int SEED = 1111;

int ElementsN = 3; // big face grid is ElementsN x ElementsN
int MIN_STEPS = 100; // minimum random walk number of steps (actual number is more to find back the start position)
float cubeL = 210; // cube size
boolean SHOW_BOX = true; // show big box frame

Element [] elements = new Element[ElementsN*ElementsN]; // small faces

// one small face (class could be used for something that doesn't look like a face)
class Element
{
  float seed = random(100,10000);
  
  ArrayList<PVector> positions = new ArrayList<PVector>();
  ArrayList<Integer> rotx = new ArrayList<Integer>();
  ArrayList<Integer> roty = new ArrayList<Integer>();
  ArrayList<Float> orientation = new ArrayList<Float>();
  
  int elemIndex;
  int nextIndex;
  
  int mParticles = 140; // number of particles on a small face
  
  Element(int k)
  {
    elemIndex = k;
  }
  
  // find element index of the destination of the face (at the end of the loop)
  void findNextElemIndex()
  {
    PVector endPos = positions.get(positions.size()-1);
    int res = 0;
    float bestDistance = 12345;
    for(int ind=0;ind<(ElementsN*ElementsN-1);ind++)
    {
      PVector startPos = elements[ind].positions.get(0);
      float d = dist(ElementsN-1-endPos.x,ElementsN-1-endPos.y,startPos.x,startPos.y);
      if(d<bestDistance)
      {
        bestDistance = d;
        res = ind;
      }
    }
    nextIndex = res;
  }
  
  PVector particleNoise(int i, float p)
  {
    float scl = 2.0;
    float nx = (float)noise.eval((i+1)*seed, scl*p);
    float ny = (float)noise.eval((i+1)*seed, 1000+scl*p);
    float nz = 2.0*pow(map((float)noise.eval((i+1)*seed, 2000+4*scl*p),-1,1,0,1),3.1); // used in size of dot
    return new PVector(nx,ny,nz);
  }
  
  // lerping towards destination noise
  PVector combinedParticleNoise(int i, float p)
  {
    PVector ns1 = particleNoise(i,p);
    Element endElem = elements[nextIndex];
    PVector ns2 = endElem.particleNoise(i,p-1); // noise of next face seen at its p-1
    float rx = rotx.get(rotx.size()-1);
    float ry = roty.get(roty.size()-1);
    float sgnx = rx>0?1:-1;
    float sgny = ry>0?1:-1;
    float sgn3 = (rx+ry)%2==0?1:-1;
    PVector ns3 = new PVector(sgn3*sgnx*ns2.x, sgn3*sgny*ns2.y,ns2.z);
    PVector ns = new PVector(lerp(ns1.x,ns3.x,p), lerp(ns1.y,ns3.y,p), lerp(ns1.z,ns3.z,p));
    return ns;
  }
    
  float pToFloatIndex(float p)
  {
    return (elements[0].positions.size()-1)*p;
  }

  // draw a face of size a
  void show0(float p, float a)
  {
    push();
    
    int index = int(floor(pToFloatIndex(p)));
    
    float rx = rotx.get(index);
    float ry = roty.get(index);
    
    rotateX(PI*rx);
    rotateY(PI*ry);
    
    for(int i=0;i<mParticles;i++)
    {
      PVector ns = combinedParticleNoise(i,p);
      
      float S = 2.5;
      float x = a/2*tanh(ns.x*S); // tanh bounding (particles stay in square)
      float y = a/2*tanh(ns.y*S);
      
      float sz = 0.45+1.5*ns.z;
      
      // size propagation on some particles
      if(i%2>0)
      {
        float offset = -0.003*modelZ(x,y,0);
        sz *= 2.2*pow(map(sin(TWO_PI*(4*t-offset)),-1,1,0.4,1),4.0);
      }
      fill(255);
      noStroke();
      sphereDetail(5);
      
      push();
      translate(x,y);
      sphere(sz);  
      pop();
    }
    
    float sphereSize = 2;
    fill(255);
    noStroke();
        
    push();
    translate(-a/2,-a/2);
    sphere(sphereSize);
    pop();
    
    push();
    translate(a/2,-a/2);
    sphere(sphereSize);
    pop();
    
    push();
    translate(a/2,a/2);
    sphere(sphereSize);
    pop();
    
    push();
    translate(-a/2,a/2);
    sphere(sphereSize);
    pop();
    
    /* debug
    float or = orientation.get(index);
    fill(or>0?color(0,255,0):color(0,0,255));
    */
    fill(0);
    stroke(200);
    strokeWeight(1.6);
    
    rectMode(CENTER);
    rect(0,0,a,a);
    
    pop();
  }
  
void show(float p)
  {
    float indf = pToFloatIndex(p);
    int i1 = floor(indf);
    int i2 = i1+1;
    float frac = indf-i1;
    
    PVector pos1 = positions.get(i1);
    PVector pos2 = positions.get(i2);
    
    PVector pos = pos1.copy().lerp(pos2,ease(frac,1.7)); // try other easing?
    
    float smallL = cubeL/ElementsN;
    
    // from integer indices to pixel positions
    pos.add(new PVector(-ElementsN/2.0+0.5,-ElementsN/2.0+0.5));
    pos.mult(smallL);
    PVector pixelPos1 = pos1.copy();
    PVector pixelPos2 = pos2.copy();
    pixelPos1.add(new PVector(-ElementsN/2.0+0.5,-ElementsN/2.0+0.5));
    pixelPos1.mult(smallL); 
    pixelPos2.add(new PVector(-ElementsN/2.0+0.5,-ElementsN/2.0+0.5));
    pixelPos2.mult(smallL);
    
    PVector dir = pixelPos2.copy().sub(pixelPos1);
    
    if(dir.mag()==0) // no flip
    {
      push();
      translate(pos.x,pos.y);
      show0(p,smallL);
      pop();
    }
    else // flip
    {
      push();
      translate(pixelPos1.x,pixelPos1.y);
      
      float easedFrac = ease(frac,2.2);
      int index = int(floor(pToFloatIndex(p)));
      float o = orientation.get(index);
      
      float rotSign1 = o>0?1:-1;
      float rotSign2 = dir.x-dir.y>0?1:-1;
      float rotSign = rotSign1*rotSign2;
      
      if(abs(dir.x)!=0)
      {
        translate(dir.x/2,0,0);
        rotateY(rotSign*easedFrac*PI);
        translate(-dir.x/2,0,0);
      }
      else
      {
        translate(0,dir.y/2,0);
        rotateX(rotSign*easedFrac*PI);
        translate(0,-dir.y/2,0);
      }
      show0(p,smallL);
      pop();
    }
  }
}

// Show large face
void showGrid(float p)
{
  for(int e=0;e<(ElementsN*ElementsN-1);e++)
  {
    elements[e].show(p);
  }
}

// easing on big faces
float mp01_2(float x,float a,float b)
{
  return ease(c01(3.35*mp01(x,a,b)),2.8);
}

void mainRotation()
{
  rotateX(-0.393*HALF_PI);
  rotateY(0.5*HALF_PI);
}

// factorization for use in both big faces drawing and small faces flip orientations search
void pushLargeStructureTransformations(float p,boolean full)
{
  if(full)
  {
    push();
    mainRotation();
  }
  
  //float q = pow(1-p,1.3); // easing is possible
  float q = p; // here no easing on global time
  
  float p1 = mp01_2(q,0./6,1.0/6);
  float p2 = mp01_2(q,1.0/6,2.0/6);
  float p3 = mp01_2(q,2.0/6,3.0/6);
  float p4 = mp01_2(q,3.0/6,4.0/6);
  float p5 = mp01_2(q,4.0/6,5.0/6);
  float p6 = mp01_2(q,5.0/6,6.0/6);
  
  push();
  translate(0,0,cubeL/2);

  push();
  translate(0,-cubeL/2,0);
  rotateX(3*p1*HALF_PI);
  translate(0,+cubeL/2,0);
  
  push();
  translate(-cubeL/2,0,0);
  rotateY(3*p2*HALF_PI);
  translate(+cubeL/2,0,0);
  
  push();
  translate(0,cubeL/2,0);
  rotateX(-3*p3*HALF_PI);
  translate(0,-cubeL/2,0);
  
  push();
  translate(cubeL/2,0,0);
  rotateY(-3*p4*HALF_PI);
  translate(-cubeL/2,0,0);
  
  push();
  translate(0,-cubeL/2,0);
  rotateX(3*p5*HALF_PI);
  translate(0,+cubeL/2,0);
  
  push();
  translate(-cubeL/2,0,0);
  rotateY(3*p6*HALF_PI);
  translate(+cubeL/2,0,0);
}

void popLargeStructureTransformations(float p,boolean full)
{
  pop();
  pop();
  pop();
  pop();
  pop();
  pop();
  pop();
  
  if(full) pop();
}

void drawStructure(float p)
{
  p = (p+1234)%1;
  
  push();
  mainRotation();
  
  if(SHOW_BOX)
  {
    stroke(200);
    strokeWeight(0.8);
    noFill();
    box(cubeL);
  }
  
  for(int a=-1;a<=1;a+=2)
  {
    for(int b=-1;b<=1;b+=2)
    {
      for(int c=-1;c<=1;c+=2)
      {
        float x = cubeL*a/2;
        float y = cubeL*b/2;
        float z = cubeL*c/2;
        push();
        translate(x,y,z);
        fill(255);
        noStroke();
        sphere(3.5);
        pop();
      }
    }
  }
  
  pushLargeStructureTransformations(p,false);
  
  showGrid(p);
  
  popLargeStructureTransformations(p,false);
  
  pop();
}

// random walk that avoids going back to previous position
void buildRandomPath()
{
  PVector holePos = new PVector(1,1);
  int [][] indexAtPos = new int[ElementsN][ElementsN];

  int k = 0;
  for(int i=0;i<ElementsN;i++)
  {
    for(int j=0;j<ElementsN;j++)
    {
      if(holePos.x!=i || holePos.y!=j)
      {
        elements[k] = new Element(k);
        elements[k].positions.add(new PVector(i,j));
        indexAtPos[i][j] = k;
        k++;
      }
    }
  }
  
  PVector startHolePos = holePos;
  PVector comingFrom = new PVector(-12345,-12345);
  
  int iteration = 0;
  
  while(iteration<=MIN_STEPS||holePos.x!=(ElementsN-1-startHolePos.x)||holePos.y!=(ElementsN-1-startHolePos.y))
  {
    boolean notok = true;
    PVector previousHolePos = holePos;

    while(notok)
    {
      int choice = floor(random(4));
      
      PVector target = new PVector(-12345,-12345);
      
      if(choice==0) target = new PVector(holePos.x-1,holePos.y);
      if(choice==1) target = new PVector(holePos.x+1,holePos.y);
      if(choice==2) target = new PVector(holePos.x,holePos.y-1);
      if(choice==3) target = new PVector(holePos.x,holePos.y+1);
      
      if(target.x==comingFrom.x&&target.y==comingFrom.y)
          continue;
      if(target.x<0||target.x>=ElementsN||target.y<0||target.y>=ElementsN)
          continue;
          
      PVector aux = holePos;
      holePos = target;
      previousHolePos = aux;
      notok = false;
    }
    
    int elementInd = indexAtPos[round(holePos.x)][round(holePos.y)];
    
    indexAtPos[round(previousHolePos.x)][round(previousHolePos.y)] = elementInd;

    elements[elementInd].positions.add(previousHolePos);
    
    for(int e=0;e<(ElementsN*ElementsN-1);e++)
    {
      if(e!=elementInd)
      {
        elements[e].positions.add(elements[e].positions.get(iteration));
      }
    }
    
    comingFrom = previousHolePos;
    
    iteration++;
  }
  
  // Looking for orientation to flip the small faces outside the box
  // also recording done rotations due to flips
  // and looking for destination index in elements array
  for(int e=0;e<(ElementsN*ElementsN-1);e++)
  {
    int index = 0;
    int curRotx = 0;
    int curRoty = 0;
    for(PVector pos : elements[e].positions)
    {
      float p = 1.0*(index+1)/elements[e].positions.size();
      elements[e].rotx.add(curRotx);
      elements[e].roty.add(curRoty);
      
      pushLargeStructureTransformations(p,true);
      
      PVector v1 = new PVector(modelX(0,0,0),modelY(0,0,0),modelZ(0,0,0)); // position and vector to big face
      PVector v0 = new PVector(modelX(0,0,10),modelY(0,0,10),modelZ(0,0,10)); // position in front of big face
      PVector v2 = v0.copy().sub(v1);
      
      v1.normalize();
      v2.normalize();
      float scalarProduct = v1.dot(v2);
      elements[e].orientation.add(scalarProduct);
      
      popLargeStructureTransformations(p,true);
      
      index++;
      
      // recording done rotations due to flips
      if(index<elements[e].positions.size()) curRotx = (curRotx+int(abs(elements[e].positions.get(index).y-elements[e].positions.get(index-1).y)))%2;
      if(index<elements[e].positions.size()) curRoty = (curRoty+int(abs(elements[e].positions.get(index).x-elements[e].positions.get(index-1).x)))%2;
    }
    
    // find index of destination in elements array
    elements[e].findNextElemIndex();
  }
}

void setup(){
  size(600,600,P3D);
  result = new int[width*height][3];
  
  randomSeed(SEED);
  
  noise = new OpenSimplexNoise();
  
  ortho();

  buildRandomPath();
}


void draw_(){
  background(0);
  push();
  translate(width/2,height/2);
  
  //scale(0.9);
  
  // replacement technique, showing 5 big faces
  int K=5;
  for(int i=0;i<K;i++)
  {
    float p = (i+t)/K;
    drawStructure(p);
  }

  pop();
}
