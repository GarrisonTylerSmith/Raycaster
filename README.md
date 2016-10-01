# Raycaster
I am creating a new project in which I will be creating a basic raycaster . I will also have to read in a json file which I will try have modified based on my project
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
int line = 1;

// next_c() wraps the getc() function and provides error checking and line
// number maintenance
int next_c(FILE* json) {
  int c = fgetc(json);
#ifdef DEBUG
  printf("next_c: '%c'\n", c);
#endif
  if (c == '\n') {
    line += 1;
  }
  if (c == EOF) {
    fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", line);
    exit(1);
  }
  return c;
}


// expect_c() checks that the next character is d.  If it is not it emits
// an error.
void expect_c(FILE* json, int d) {
  int c = next_c(json);
  if (c == d) return;
  fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, line);
  exit(1);    
}


// skip_ws() skips white space in the file.
void skip_ws(FILE* json) {
  int c = next_c(json);
  while (isspace(c)) {
    c = next_c(json);
  }
  ungetc(c, json);
}


// next_string() gets the next string from the file handle and emits an error
// if a string can not be obtained.
char* next_string(FILE* json) {
  char buffer[129];
  int c = next_c(json);
  if (c != '"') {
    fprintf(stderr, "Error: Expected string on line %d.\n", line);
    exit(1);
  }  
  c = next_c(json);
  int i = 0;
  while (c != '"') {
    if (i >= 128) {
      fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
      exit(1);      
    }
    if (c == '\\') {
      fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
      exit(1);      
    }
    if (c < 32 || c > 126) {
      fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
      exit(1);
    }
    buffer[i] = c;
    i += 1;
    c = next_c(json);
  }
  buffer[i] = 0;
  return strdup(buffer);
}

double next_number(FILE* json) {
  double value;
  fscanf(json, "%f", &value);
  // Error check this..
  return value;
}

double* next_vector(FILE* json) {
  double* v = malloc(3*sizeof(double));
  expect_c(json, '[');
  skip_ws(json);
  v[0] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[1] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[2] = next_number(json);
  skip_ws(json);
  expect_c(json, ']');
  return v;
}


void read_scene(char* filename) {
  int c;
  FILE* json = fopen(filename, "r");

  if (json == NULL) {
    fprintf(stderr, "Error: Could not open file \"%s\"\n", filename);
    exit(1);
  }
  
  skip_ws(json);
  
  // Find the beginning of the list
  expect_c(json, '[');

  skip_ws(json);

  // Find the objects

  while (1) {
    c = fgetc(json);
    if (c == ']') {
      fprintf(stderr, "Error: This is the worst scene file EVER.\n");
      fclose(json);
      return;
    }
    if (c == '{') {
      skip_ws(json);
    
      // Parse the object
      char* key = next_string(json);
      if (strcmp(key, "type") != 0) {
	fprintf(stderr, "Error: Expected \"type\" key on line number %d.\n", line);
	exit(1);
      }

      skip_ws(json);

      expect_c(json, ':');

      skip_ws(json);

      char* value = next_string(json);

      if (strcmp(value, "camera") == 0) {
      } else if (strcmp(value, "sphere") == 0) {
      } else if (strcmp(value, "plane") == 0) {
      } else {
	fprintf(stderr, "Error: Unknown type, \"%s\", on line number %d.\n", value, line);
	exit(1);
      }

      skip_ws(json);

      while (1) {
	// , }
	c = next_c(json);
	if (c == '}') {
	  // stop parsing this object
	  break;
	} else if (c == ',') {
	  // read another field
	  skip_ws(json);
	  char* key = next_string(json);
	  skip_ws(json);
	  expect_c(json, ':');
	  skip_ws(json);
	  if ((strcmp(key, "width") == 0) ||
	      (strcmp(key, "height") == 0) ||
	      (strcmp(key, "radius") == 0)) {
	    double value = next_number(json);
	  } else if ((strcmp(key, "color") == 0) ||
		     (strcmp(key, "position") == 0) ||
		     (strcmp(key, "normal") == 0)) {
	    double* value = next_vector(json);
	  } else {
	    fprintf(stderr, "Error: Unknown property, \"%s\", on line %d.\n",
		    key, line);
	    //char* value = next_string(json);
	  }
	  skip_ws(json);
	} else {
	  fprintf(stderr, "Error: Unexpected value on line %d\n", line);
	  exit(1);
	}
      }
      skip_ws(json);
      c = next_c(json);
      if (c == ',') {
	// noop
	skip_ws(json);
      } else if (c == ']') {
	fclose(json);
	return;
      } else {
	fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", line);
	exit(1);
      }
    }
  }
}

int main(int c, char** argv) {
  read_scene(argv[1]);
  return 0;
}
The Code above is the Parsing Json file which allows me to parse the Json file . The only thing I need to change on it is to create a if statement to make sure that I am checking that there is not a radius for a plane.



// polymorphism in c
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
typedef struct{
	int kind; // 0 = cylinder, 1 = sphere, 2 = teapot
	double color[3];
	union{
		struct {
			double center[3];
			double radius;
			double height;
		}cylinder;
		struct{
			double center[3];
			double radius;
		}sphere;
		struct{
			double handle_length;
		}teapot;
		struct{
			double position;
			double normal;
		}plane;
		struct{
			double position;
			double width;
			double height;
		}camera;
	};
}Object;

static inline double sqr(double v){
	return v*v;
}
static inline void normalize(double* v){
	double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
	v[0] /= len;
	v[1] /= len;
	v[2] /= len;
}
double cylinder_intersecion(double* Ro, double* Rd, double* c, double r){
	// Step 1 . Find the equation for the object 
	// we are interested in (this case a cylinder)
	// this is essential for the camera view point
	// x^2 +z^2 = r^2
	// Step 2. Parameterize the equation with a center point (if needeed)
	// (x-Cx)^2 + (z-Cz)^2 = r^2
	// step 3. Substitute the eq for a ray into our object eq.
	//
	//(RoX + RDX*t - Cx)^2 + ( Roz + t*Rdz - Cz)^2 - r^2 = 0
	//
	// Step 4. Solve for t .
	//
	//Step 4.a Rewrite the equation(flatten).
	//
	// -r^2 + t^2 * Rdx^2 + t^2 * Rdz^2 + 2*t * Rox * Rdx - 2*t * Rdx * Cx + 2*t * Roz * Rdz - 2*t * Rdz *cz + Rox^2 - 2*Rox*cx + cx^2 + Roz^2 - 2* Roz* Cz + Cz^2 = 0
	//
	// Step 4.b Rewrite the equation in terms of t.
	//
	// t^2 * (Rdx^2 + Rdz^2) + 
	// t* (2 * (Rox * Rdx - Rdx * cx + Roz * Rdz - Rdz * Cz)) + 
	// Rox^2 - 2*Rox*cx + cx^2 + Roz^2 - 2* Roz* Cz + Cz^2 - r^2 = 0
	// Use quadratic equation to solvew for t..
	//
	double a = (sqr(Rd[0]) + sqr(Rd[2]));
	double b = (2 * (Ro[0] * Rd[0] - Rd[0] * c[0] + Ro[2] * Rd[2] - Rd[2] * c[2]))
	double c = (sqr(Ro[0]) - 2*Ro[0]*c[0] + sqr(c[0]) + sqr(Ro[2]) - 2* Ro[2]*c[2] + sqr(c[2]) - sqr(r));
	double det = sqr(b) - 4*a*c; 
	if( det < 0 ) return -1;

	det = sqrt(det);

	double t0 = (-b - det) / (2*a);
	if (t0 > 0 ) return t0;

	double t1 = (-b + det) / (2*a);
	if(t1 > 0) return t1;

	return -1;

}
double ray_plane_intersection(double* Ro, double* Rd, double D, double* n){
	return -(n[0] + n[1] + n[2] + D) / (Rd[0] + Rd[1] + Rd[2]);
	//Ray: R(t) = Ro + Rd*t
	//Plane: ax + by + cz + d = 0
	//solve for t
	//n is the normal plane
	// D is the distance from the origin
	// do subsituting we get 
	// t = -(n*Ro + d)/ (n * Rd) 

	// double a = -(n*Ro + d)
	// double b = (n*Rd)
	// double det = (a/b)
	// if(det < 0) return -1;
}
double ray_sphere_intersection(double* Ro, double* Rd, double* c , double r,double* t){
	// Ray: P = Ro + Rd*t
	// Sphere: (x-xc)^2 + (y-yc)^2 + (z-zc)^2 - r^2 = 0
	// Subsituting R(t) into sphere equation
	// (Xo + Xd*t - Xc)^2 + (Yo + Yd*t - Yc)^2 + (Zo + Zd*t - Zc)^2 - r^2 = 0
	// rearranging we get
	// At^2 + Bt + C = 0 where we get two solutions
	// this is a vector from p to c
	double* vpc = c-Ro; 
	// when the sphere is behind the origin p 
	// note that this case may be dismissed if it is
	// considered that p is outside the sphere
	if((vpc * d) < 0)
		if(|vpc| > r){
			// there is no intersection
		} else if(|vpc| == r){
			intersection = p;
		}else{
			// occurs when p is inside the sphere
			pc = projecti
		}
	return -1;
}
void raycaster(double* Ro, double* Rd){
	double cx;
	double cy;
	double h;
	double w;
	// width in pixels
	double N;
	// height in pixels 
	double M;
	Pix_height = h/M;
	Pix_width = w/N;
	P_z = -Zp;
	for(int i = 0; i < M; i+=1){
		for(int j = 0; j < N; j+=1){
			double Ro[3] = {0,0,0}
			double Rd[3]{
				P_y = cy - h/2.0 + Pix_height * (i + .5),
				P_y = cx - w/2.0 + Pix_width * (j + .5),
				1
			};
			normalize(Rd); 
		}
	}

	// here I must create my raycaster
}
// create my own raycasting fucntion
int main(){
// this is all hard code, use the parsing file for json and raycaster function
	Object** objects;
	objects = malloc(sizeof(Object*)*2);
	objects[0] = malloc(sizeof(Object));
	objects[0]->kind = 0;
	objects[0]->cylinder.radius = 2;
	objects[0]->cylinder.center[0] = 0;
	objects[0]->cylinder.center[1] = 0;
	objects[0]->cylinder.center[2] = 20;
	objects[1] = NULL;

	double cx = 0;
	double cy = 0;
	double h = 0.7;
	double w = 0.7;

	int M = 20;
	int N = 20;


	double pixheight = h / M;
	double pixwidth = w / N;
	for(int y = 0; y < M; y+=1){
		for(int x = 0; x < N; x+=1){
			double Ro[3] = {0,0,0};
			// Rd = normalize(P - Ro)
			double Rd[3] = {
				cx - (w/2) + pixwidth * (x + 0.5),
				cy - (h/2) + pixheight * (y + 0.5),
				1
			};
			normalize(Rd);

			double best_t = INFINITY;
			for(int i = 0; objects[i] != 0; i+=1){
				double t = 0;
				switch(objects[i]->kind){
				case 0;
					t = cylinder_intersecion(Ro,Rd,objects[i]->cylinder.center, objects[i]->cylinder.radius);


					break;
				default;
				// horrible error
				exit(1);
				}
				if(t > 0 best_t && t < best_t) best_t = t;
			}
			if(best_t > 0 && best_t != INFINITY){
				printf("#");
			}else {
				printf(".");
			}
		}
	}


	return 0;
}

So Far in this code it is the ray tracer which I have done the plane and I am working on the sphere I had kept the cylinder which we had worked on in class because it sounds cool to put that in the view point. I am still working on the next to functions that would complete this which are the raycster function and the sphere intersection function. 
