typedef float data_t;

typedef struct {
	scalar_t x;
	scalar_t y;
} vector_t;

typedef struct {
	scalar_t density;
	scalar_t length;
	scalar_t radius;
	scalar_t elastic_modulus;
	vector_t end[2];
} rod_t;

typedef struct {
	rod_t * elements;
	int n_elements;
	scalar_t length;
} tine_t;

typedef struct {
	
	
	
} support_t;



tine_t * tine_create(scalar_t length, int n_elements, scalar_t density, scalar_t radius, int taper, int tip_weight)
{
	int k;
	time_t * tine;
	
	tine = malloc(sizeof(tine_t));
	
	tine->elements = malloc(n_elements * sizeof(rod_t));
	
	for (k = 0; k < n_elements; k++){
		tine->elements[k].density = density;
		tine->elements[k].length = length / (scalar_t) n_elements;
		tine->elements[k].radius = 
		
	}
	
	
}


	


	