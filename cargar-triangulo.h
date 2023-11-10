
typedef struct punto
{
float x, y, z, u,v;
} punto;

typedef struct vector
{
float x, y, z;
} vector;

typedef struct hiruki
{
punto p1,p2,p3;
vector v_normal;
} hiruki;

int cargar_triangulos(char *fitxiz, int *hkopptr, hiruki **hptrptr);
int cargar_triangulos_color(char *fitxiz, int *hkopptr, hiruki **hptrptr, unsigned char **rgbptr);
