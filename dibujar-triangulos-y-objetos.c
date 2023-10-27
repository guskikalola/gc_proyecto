//	Program developed by
//
//	Informatika Fakultatea
//	Euskal Herriko Unibertsitatea
//	http://www.ehu.eus/if
//
// to compile it: gcc dibujar-triangulos-y-objetos.c -lGL -lGLU -lglut
//
//
//

#include <GL/glut.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cargar-triangulo.h"

#define DESPLAZAMIENTO_TRANSLACION 5
#define ANGULO_ROTACION 3.14159 / 90
#define PROPORCION_ESCALADO 1.2

#define TRANSLACION 't'
#define ESCALADO 's'
#define ROTACION 'r'

#define SISTEMA_MUNDO 0
#define SISTEMA_LOCAL 1

#define EJE_X 0
#define EJE_Y 1
#define EJE_Z 2
#define EJE_NULL 3

#define DIR_ADELANTE 0
#define DIR_ATRAS 1

typedef struct mlist
{
    double m[16];
    struct mlist *hptr;
} mlist;

typedef struct triobj
{
    hiruki *triptr;
    int num_triangles;
    mlist *mptr;
    struct triobj *hptr;
} triobj;

// testuraren informazioa
// información de textura

extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int *dimyptr);
unsigned char *bufferra;
int dimx, dimy;

int indexx;
hiruki *triangulosptr;
triobj *foptr;
triobj *sel_ptr;
triobj *camara_ptr; // Puntero a la camara activa
int camara_activa;  // Si tenemos la camara activa o no
int denak;
int lineak;
int objektuak;
char aldaketa;
int ald_lokala;

char fitxiz[100];

// TODO
// funtzio honek u eta v koordenatuei dagokien pointerra itzuli behar du.
// debe devolver el pointer correspondiente a las coordenadas u y v
unsigned char *color_textura(float u, float v)
{
    char *lag;

    int xp, yp;

    xp = u * dimx;
    yp = (1 - v) * dimy; // Hay que invertir la v porque la y va al reves ( en el array crece de arriba a abajo; en la imagen de abajo a arriba )

    lag = (unsigned char *)bufferra;       // pixel on the left and top
    return (lag + dimx * yp * 3 + xp * 3); // Hay que multiplicar por 3 porque cada punto son 3 posiciones (r,g,b)
}

void dibujar_linea(punto p1, punto p2)
{
    float j;
    unsigned char r, g, b;
    unsigned char *colorv;
    punto *pcortemayor, *pcortemenor;
    punto pcalculado;
    double cambioj;
    if (p1.x > p2.x)
    {
        pcortemayor = &p1;
        pcortemenor = &p2;
    }
    else
    {
        pcortemayor = &p2;
        pcortemenor = &p1;
    }

    if (pcortemayor->x - pcortemenor->x == 0)
        cambioj = 1;
    else
        cambioj = 1 / (pcortemayor->x - pcortemenor->x);

    for (j = 1; j > 0; j -= cambioj)
    {
        pcalculado.x = j * pcortemayor->x + (1 - j) * pcortemenor->x;
        pcalculado.y = j * pcortemayor->y + (1 - j) * pcortemenor->y;
        pcalculado.z = j * pcortemayor->z + (1 - j) * pcortemenor->z;
        pcalculado.u = j * pcortemayor->u + (1 - j) * pcortemenor->u;
        pcalculado.v = j * pcortemayor->v + (1 - j) * pcortemenor->v;

        glBegin(GL_POINTS);
        colorv = color_textura(pcalculado.u, pcalculado.v);
        r = colorv[0];
        g = colorv[1];
        b = colorv[2];
        glColor3ub(r, g, b);
        glVertex3f(pcalculado.x, pcalculado.y, pcalculado.z);
        glEnd();
    }
}

void print_matrizea(char *str)
{
    int i;

    printf("%s\n", str);
    for (i = 0; i < 4; i++)
        printf("%lf, %lf, %lf, %lf\n", sel_ptr->mptr->m[i * 4], sel_ptr->mptr->m[i * 4 + 1], sel_ptr->mptr->m[i * 4 + 2],
               sel_ptr->mptr->m[i * 4 + 3]);
}

// TODO
// aurrerago egitekoa
// para más adelante
void mxp(punto *pptr, double m[16], punto p)
{
    // print_matrizea("m:");

    pptr->x = m[0] * p.x + m[1] * p.y + m[2] * p.z + m[3];
    pptr->y = m[4] * p.x + m[5] * p.y + m[6] * p.z + m[7];
    pptr->z = m[8] * p.x + m[9] * p.y + m[10] * p.z + m[11];
    pptr->u = p.u;
    pptr->v = p.v;
}

void mxm(double mresptr[16], double mA[16], double mB[16])
{
    /*
    [ 0  1  2  3  |       [ 0  1  2  3  |
    | 4  5  6  7  |  \/   | 4  5  6  7  |
    | 8  9  10 11 |  /\   | 8  9  10 11 |
    | 12 13 14 15 ]       | 12 13 14 15 ]

    */

    int iA, iB;
    for (iA = 0; iA <= 12; iA += 4)
    {
        for (iB = 0; iB <= 3; iB++)
        {
            mresptr[iA + iB] = mA[iA] * mB[iB] + mA[iA + 1] * mB[iB + 4] + mA[iA + 2] * mB[iB + 8] + mA[iA + 3] * mB[iB + 12];
        }
    }
}

void aplicar_mcsr(double mresptr[16], double mCamara[16], double mObj[16])
{
    double mcsr[16];
    int i;
    for (i = 0; i < 16; i++)
        mcsr[i] = 0;

    /*
        m[16] posiciones
        [ 0  1  2  3  |
        | 4  5  6  7  |
        | 8  9  10 11 |
        | 12 13 14 15 ]
    */

    /*
        Mcsr =
        [ Xc_x  Xc_y  Xc_z  (- Xc * Pc)  |
        | Yc_x  Yc_y  Yc_z  (- Yc * Pc)  |
        | Zc_x  Zc_y  Zc_z  (- Zc * Pc)  |
        |   0     0     0          1     ]

        Mcamara =
        [ Xc_x  Yc_x  Zc_x  Pc_x  |
        | Xc_y  Yc_y  Zc_y  Pc_y  |
        | Xc_z  Yc_z  Zc_z  Pc_z  |
        |   0     0     0    1    ]
    */
    // Columna Xc
    mcsr[0] = mCamara[0];
    mcsr[1] = mCamara[4];
    mcsr[2] = mCamara[8];

    // Columna Yc
    mcsr[4] = mCamara[1];
    mcsr[5] = mCamara[5];
    mcsr[6] = mCamara[9];

    // Columna Zc
    mcsr[8] = mCamara[2];
    mcsr[9] = mCamara[6];
    mcsr[10] = mCamara[10];

    // Columna -(X,Y,Z)c * Pc
    mcsr[3] = -(mCamara[0] * mCamara[3] + mCamara[4] * mCamara[7] + mCamara[8] * mCamara[11]);   // -Xc * Pc
    mcsr[7] = -(mCamara[1] * mCamara[3] + mCamara[5] * mCamara[7] + mCamara[9] * mCamara[11]);   // -Yc * Pc
    mcsr[11] = -(mCamara[2] * mCamara[3] + mCamara[6] * mCamara[7] + mCamara[10] * mCamara[11]); // -Zc * Pc

    mcsr[15] = 1;

    mxm(mresptr, mcsr, mObj);
}

void dibujar_triangulo(triobj *optr, int i)
{
    hiruki *tptr;

    punto *pgoiptr, *pbeheptr, *perdiptr;
    punto p1, p2, p3;

    float t, s;
    float cambiot, cambios;

    punto pcorte1, pcorte2;

    if (i >= optr->num_triangles)
        return;
    tptr = optr->triptr + i;

    // Si tenemos camara multiplicar la 'Mcsr' a la matriz del objeto 'Mobj'
    mlist *matriz_observador_ptr;
    if (camara_ptr != 0 && camara_activa == 1)
    {
        mlist mcsr_objptr;
        matriz_observador_ptr = &mcsr_objptr;
        aplicar_mcsr(matriz_observador_ptr->m, camara_ptr->mptr->m, optr->mptr->m);
    }
    else
    {
        matriz_observador_ptr = optr->mptr;
    }

    // (Mcsr * Mobj) * Obj
    mxp(&p1, matriz_observador_ptr->m, tptr->p1);
    mxp(&p2, matriz_observador_ptr->m, tptr->p2);
    mxp(&p3, matriz_observador_ptr->m, tptr->p3);

    if (lineak == 1)
    {
        glBegin(GL_POLYGON);
        glVertex3d(p1.x, p1.y, p1.z);
        glVertex3d(p2.x, p2.y, p2.z);
        glVertex3d(p3.x, p3.y, p3.z);
        glEnd();
        return;
    }
    //  else
    // TODO
    // hemen azpikoa kendu eta triangelua testurarekin marrazten duen kodea sartu.0
    // lo que sigue aqui hay que sustituir por el código adecuado que dibuja el triangulo con textura

    // Calcular Psup, Pmed, Pinf

    if (p1.y > p2.y)
    {
        pgoiptr = &p1;  // Psup <- P1
        pbeheptr = &p2; // Pinf <- P2
    }
    else
    {
        pgoiptr = &p2;  // Psup <- P2
        pbeheptr = &p1; // Pinf <- P1
    }

    if (p3.y > pgoiptr->y)
    {
        perdiptr = pgoiptr; // Pmed <- Psup
        pgoiptr = &p3;      // Psup <- P3
    }
    else if (p3.y < pbeheptr->y)
    {
        perdiptr = pbeheptr; // Pmed <- Pinf
        pbeheptr = &p3;      // Psup <- P3
    }
    else
        perdiptr = &p3; // Pmed <- P3

    // Triangulo es una linea
    if (p1.y == p2.y && p2.y == p3.y)
    {
        // Coger la y min y la y max y punto de min a max
        if (p1.x > p2.x)
        {
            pcorte1 = p1;
            pcorte2 = p2;
        }
        else
        {
            pcorte1 = p2;
            pcorte2 = p1;
        }

        if (p3.x > pcorte1.x)
        {
            pcorte1 = p3;
        }
        else if (p3.x < pcorte2.x)
        {
            pcorte2 = p3;
        }

        dibujar_linea(pcorte1, pcorte2);
        return;
    }

    if (pgoiptr->y - perdiptr->y == 0)
        cambiot = 1;
    else
        cambiot = 1 / (pgoiptr->y - perdiptr->y);

    if (pgoiptr->y - pbeheptr->y == 0)
        cambios = 1;
    else
        cambios = 1 / (pgoiptr->y - pbeheptr->y);

    // Vamos a utilizar las coordenadas baricentricas para calcular (x,y,z,u,v) para
    // cada punto interior del triangulo. Para ello vamos a dividir el triangulo en
    // dos secciones, la superior y la inferior. De esta manera solo tenemos que
    // preocuparnos de ir calculando el punto de corte de dos lineas.

    // De A -> B (avanzamos con t)  y  A -> C (avanzamos con s)
    for (t = 1, s = 1; t > 0; t = t - cambiot, s = s - cambios)
    {
        pcorte1.x = t * pgoiptr->x + (1 - t) * perdiptr->x;
        pcorte1.y = t * pgoiptr->y + (1 - t) * perdiptr->y;
        pcorte1.z = t * pgoiptr->z + (1 - t) * perdiptr->z;
        pcorte1.u = t * pgoiptr->u + (1 - t) * perdiptr->u;
        pcorte1.v = t * pgoiptr->v + (1 - t) * perdiptr->v;

        pcorte2.x = s * pgoiptr->x + (1 - s) * pbeheptr->x;
        pcorte2.y = s * pgoiptr->y + (1 - s) * pbeheptr->y;
        pcorte2.z = s * pgoiptr->z + (1 - s) * pbeheptr->z;
        pcorte2.u = s * pgoiptr->u + (1 - s) * pbeheptr->u;
        pcorte2.v = s * pgoiptr->v + (1 - s) * pbeheptr->v;

        dibujar_linea(pcorte1, pcorte2);
    }

    // Tenemos que sumarle para cancelar el cambio de mas que ha hecho en la ultima iteracion.
    // Si no hacemos esto la s va a tomar valores negativos en el siguiente for porque empieza
    // un desplazamiento mas abajo del que deberia
    s += cambios;

    // Dibujar parte inferior
    // Seguimos de la s anterior, es decir
    // seguimos desde el final del segmento
    // superior del triangulo.
    // La t la tenemos que volver a asignar
    // a 1 porque la linea correspondiente a t
    // la hemos cambiado.

    // De B -> C (avanzamos con t)  y  A -> C (avanzamos con s)

    if (perdiptr->y - pbeheptr->y == 0)
        cambiot = 1;
    else
        cambiot = 1 / (perdiptr->y - pbeheptr->y);

    for (t = 1; t > 0; t = t - cambiot, s = s - cambios)
    {
        pcorte1.x = t * perdiptr->x + (1 - t) * pbeheptr->x;
        pcorte1.y = t * perdiptr->y + (1 - t) * pbeheptr->y;
        pcorte1.z = t * perdiptr->z + (1 - t) * pbeheptr->z;
        pcorte1.u = t * perdiptr->u + (1 - t) * pbeheptr->u;
        pcorte1.v = t * perdiptr->v + (1 - t) * pbeheptr->v;

        pcorte2.x = s * pgoiptr->x + (1 - s) * pbeheptr->x;
        pcorte2.y = s * pgoiptr->y + (1 - s) * pbeheptr->y;
        pcorte2.z = s * pgoiptr->z + (1 - s) * pbeheptr->z;
        pcorte2.u = s * pgoiptr->u + (1 - s) * pbeheptr->u;
        pcorte2.v = s * pgoiptr->v + (1 - s) * pbeheptr->v;

        dibujar_linea(pcorte1, pcorte2);
    }
}

static void marraztu(void)
{
    float u, v;
    int i, j;
    triobj *auxptr;
    /*
    unsigned char* colorv;
    unsigned char r,g,b;
    */

    // marrazteko objektuak behar dira
    // no se puede dibujar sin objetos
    if (foptr == 0)
        return;

    // clear viewport...
    if (objektuak == 1)
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    else
    {
        if (denak == 0)
            glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    }

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-500.0, 500.0, -500.0, 500.0, -500.0, 500.0);

    triangulosptr = sel_ptr->triptr;
    if (objektuak == 1)
    {
        if (denak == 1)
        {
            for (auxptr = foptr; auxptr != 0; auxptr = auxptr->hptr)
            {
                for (i = 0; i < auxptr->num_triangles; i++)
                {
                    dibujar_triangulo(auxptr, i);
                }
            }
        }
        else
        {
            for (i = 0; i < sel_ptr->num_triangles; i++)
            {
                dibujar_triangulo(sel_ptr, i);
            }
        }
    }
    else
    {
        dibujar_triangulo(sel_ptr, indexx);
    }
    glFlush();
}

void read_from_file(char *fitx)
{
    int i, retval;
    triobj *optr;

    // printf("%s fitxategitik datuak hartzera\n",fitx);
    optr = (triobj *)malloc(sizeof(triobj));
    retval = cargar_triangulos(fitx, &(optr->num_triangles), &(optr->triptr));
    if (retval != 1)
    {
        printf("%s fitxategitik datuak hartzerakoan arazoak izan ditut\n    Problemas al leer\n", fitxiz);
        free(optr);
    }
    else
    {
        triangulosptr = optr->triptr;
        // printf("objektuaren matrizea...\n");
        optr->mptr = (mlist *)malloc(sizeof(mlist));
        for (i = 0; i < 16; i++)
            optr->mptr->m[i] = 0;
        // Inicializa la matriz de identidad
        optr->mptr->m[0] = 1.0;
        optr->mptr->m[5] = 1.0;
        optr->mptr->m[10] = 1.0;
        optr->mptr->m[15] = 1.0;
        optr->mptr->hptr = 0; // La siguiente matriz ez "null"
        // printf("objektu zerrendara doa informazioa...\n");
        optr->hptr = foptr; // el siguiente al objeto es el que antes era el primero
        foptr = optr;       // foptr apunta al primer objeto ( el ultimo cargado )
        sel_ptr = optr;
    }
    printf("datuak irakurrita\nLecura finalizada\n");
}

void translacion(mlist *matriz_transformacion, int eje, int dir, int cantidad)
{
    int i, mul_dir, n, m, o;
    for (i = 0; i < 16; i++)
        matriz_transformacion->m[i] = 0;

    if (dir == DIR_ADELANTE)
        mul_dir = 1;
    else // dir == DIR_ATRAS
        mul_dir = -1;

    if (eje == EJE_X)
    {
        m = cantidad * mul_dir;
        n = 0;
        o = 0;
    }
    else if (eje == EJE_Y)
    {
        m = 0;
        n = cantidad * mul_dir;
        o = 0;
    }
    else // EJE_Z
    {
        m = 0;
        n = 0;
        o = cantidad * mul_dir;
    }

    matriz_transformacion->m[0] = 1;
    matriz_transformacion->m[5] = 1;
    matriz_transformacion->m[10] = 1;
    matriz_transformacion->m[15] = 1;

    matriz_transformacion->m[3] = m;  // x
    matriz_transformacion->m[7] = n;  // y
    matriz_transformacion->m[11] = o; // z
}

void rotacion(mlist *matriz_transformacion, int eje, int dir, double angulo)
{
    int i;
    float angulo_x_dir;
    for (i = 0; i < 16; i++)
        matriz_transformacion->m[i] = 0;

    if (dir == DIR_ADELANTE)
        angulo_x_dir = angulo;
    else // dir == DIR_ATRAS
        angulo_x_dir = (2 * 3.14159) - angulo;

    if (eje == EJE_X)
    {
        matriz_transformacion->m[0] = 1;
        matriz_transformacion->m[5] = cos(angulo_x_dir);
        matriz_transformacion->m[6] = -sin(angulo_x_dir);
        matriz_transformacion->m[9] = sin(angulo_x_dir);
        matriz_transformacion->m[10] = cos(angulo_x_dir);
    }
    else if (eje == EJE_Y)
    {
        matriz_transformacion->m[5] = 1;
        matriz_transformacion->m[0] = cos(angulo_x_dir);
        matriz_transformacion->m[2] = sin(angulo_x_dir);
        matriz_transformacion->m[8] = -sin(angulo_x_dir);
        matriz_transformacion->m[10] = cos(angulo_x_dir);
    }
    else // EJE_Z
    {
        matriz_transformacion->m[0] = cos(angulo_x_dir);
        matriz_transformacion->m[1] = -sin(angulo_x_dir);
        matriz_transformacion->m[4] = sin(angulo_x_dir);
        matriz_transformacion->m[5] = cos(angulo_x_dir);
        matriz_transformacion->m[10] = 1;
    }

    matriz_transformacion->m[15] = 1;
}

void escalado(mlist *matriz_transformacion, int dir, float proporcion)
{
    if (proporcion <= 0)
        return; // Directamente evitamos que se pueda pasar 0, evitando division /0

    int i;
    for (i = 0; i < 16; i++)
    {
        matriz_transformacion->m[i] = 0;
    }

    if (dir == DIR_ADELANTE)
    {
        matriz_transformacion->m[0] = proporcion;  // x
        matriz_transformacion->m[5] = proporcion;  // y
        matriz_transformacion->m[10] = proporcion; // z
        matriz_transformacion->m[15] = 1;
    }
    else // dir == DIR_ATRAS
    {
        matriz_transformacion->m[0] = 1 / proporcion;  // x
        matriz_transformacion->m[5] = 1 / proporcion;  // y
        matriz_transformacion->m[10] = 1 / proporcion; // z
        matriz_transformacion->m[15] = 1;
    }
}

void aplicar_transformacion(mlist *matriz_transformacionptr, int sistema_referencia)
{
    mlist *nueva_matrizptr = (mlist *)malloc(sizeof(mlist));

    if (sistema_referencia == SISTEMA_LOCAL)
    {
        // Multiplicar por la derecha
        mxm(nueva_matrizptr->m, sel_ptr->mptr->m, matriz_transformacionptr->m);
    }
    else // SISTEMA_MUNDO
    {
        // Multiplicar por la izquierda
        mxm(nueva_matrizptr->m, matriz_transformacionptr->m, sel_ptr->mptr->m);
    }

    nueva_matrizptr->hptr = sel_ptr->mptr;
    sel_ptr->mptr = nueva_matrizptr;
}

void tratar_transformacion(int eje, int dir)
{
    mlist matriz_transformacion;

    switch (aldaketa)
    {
    case TRANSLACION:
        translacion(&matriz_transformacion, eje, dir, DESPLAZAMIENTO_TRANSLACION);
        break;
    case ROTACION:
        rotacion(&matriz_transformacion, eje, dir, ANGULO_ROTACION);
        break;
    case ESCALADO:
        if (eje != EJE_NULL)
            return;
        escalado(&matriz_transformacion, dir, PROPORCION_ESCALADO);
        break;
    default:
        break;
    }
    aplicar_transformacion(&matriz_transformacion, ald_lokala);
}

void x_aldaketa(int dir)
{
    tratar_transformacion(EJE_X, dir);
}

void y_aldaketa(int dir)
{
    tratar_transformacion(EJE_Y, dir);
}

void z_aldaketa(int dir)
{
    tratar_transformacion(EJE_Z, dir);
}

void undo()
{
    /*

    La matriz de identidad no la podemos borrar. Esa siempre va a estar en la lista.
    Sabemos que hemos llegado a la identidad porque su hptr ( siguiente matriz ) tiene valor 0.

    mlist
        m2 -> m1 -> I -> 0

    */

    mlist *matriz_a_borrarptr;

    if (sel_ptr->mptr->hptr != 0)
    {
        matriz_a_borrarptr = sel_ptr->mptr;
        sel_ptr->mptr = sel_ptr->mptr->hptr;
        free(matriz_a_borrarptr);
    }
}

// This function will be called whenever the user pushes one key
static void teklatua(unsigned char key, int x, int y)
{
    int retval;
    int i;
    FILE *obj_file;

    switch (key)
    {
    case 13:
        if (foptr != 0) // objekturik ez badago ezer ez du egin behar
                        // si no hay objeto que no haga nada
        {
            indexx++; // azkena bada lehenengoa bihurtu
                      // pero si es el último? hay que controlarlo!
            if (indexx == sel_ptr->num_triangles)
            {
                indexx = 0;
                if ((denak == 1) && (objektuak == 0))
                {
                    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
                    glFlush();
                }
            }
        }
        break;
    case 'd':
        if (denak == 1)
            denak = 0;
        else
            denak = 1;
        break;
    case 'o':
        if (objektuak == 1)
            objektuak = 0;
        else
            objektuak = 1;
        break;
    case 'l':
        if (lineak == 1)
            lineak = 0;
        else
            lineak = 1;
        break;
    case 't':
        aldaketa = 't';
        break;
    case 'r':
        aldaketa = 'r';
        break;
    case 's':
        aldaketa = 's';
        tratar_transformacion(EJE_NULL, DIR_ADELANTE);
        break;
    case 'S':
        aldaketa = 's';
        tratar_transformacion(EJE_NULL, DIR_ATRAS);
        break;
    case 'g':
        if (ald_lokala == 1)
            ald_lokala = 0;
        else
            ald_lokala = 1;

        if (ald_lokala == SISTEMA_LOCAL)
            printf("Cambiado a sistema local\n");
        else
            printf("Cambiado al sistema del mundo\n");
        break;
    case 'x':
        x_aldaketa(DIR_ADELANTE);
        break;
    case 'y':
        y_aldaketa(DIR_ADELANTE);
        break;
    case 'z':
        z_aldaketa(DIR_ADELANTE);
        break;
    case 'X':
        x_aldaketa(DIR_ATRAS);
        break;
    case 'Y':
        y_aldaketa(DIR_ATRAS);
        break;
    case 'Z':
        z_aldaketa(DIR_ATRAS);
        break;
    case 'u':
        undo();
        break;
    case 'j':
        printf("Camara activa %d\n", camara_activa);
        if (camara_activa == 1)
            camara_activa = 0;
        else
            camara_activa = 1;
        break;
    case 'f':
        /*Ask for file*/
        printf("idatzi fitxategi izena\n");
        scanf("%s", &(fitxiz[0]));
        read_from_file(fitxiz);
        indexx = 0;
        break;
        /* case 'S':  // save to file
             printf("idatzi fitxategi izena\n");
             scanf("%s", &(fitxiz[0]));
                 if ((obj_file = fopen(fitxiz, "w")) == NULL)
                          {
                          printf("ezin fitxategia ireki\n");
                          }
                      else
                          {
                          for (i =0; i < sel_ptr->num_triangles; i++)
                             {
                             fprintf(obj_file,"t %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                                  sel_ptr->triptr[i].p1.x-250, sel_ptr->triptr[i].p1.y-250, sel_ptr->triptr[i].p1.z,
                                  sel_ptr->triptr[i].p1.u, sel_ptr->triptr[i].p1.v,
                                  sel_ptr->triptr[i].p2.x-250, sel_ptr->triptr[i].p2.y-250, sel_ptr->triptr[i].p2.z,
                                  sel_ptr->triptr[i].p2.u, sel_ptr->triptr[i].p2.v,
                                  sel_ptr->triptr[i].p3.x-250, sel_ptr->triptr[i].p3.y-250, sel_ptr->triptr[i].p3.z,
                                  sel_ptr->triptr[i].p3.u, sel_ptr->triptr[i].p3.v );
                             }
                          fclose(obj_file);
                          }
                 break; */
    case 9:             /* <TAB> */
        if (foptr != 0) // objekturik gabe ez du ezer egin behar
                        // si no hay objeto no hace nada
        {
            sel_ptr = sel_ptr->hptr;
            /*The selection is circular, thus if we move out of the list we go back to the first element*/
            if (sel_ptr == 0)
                sel_ptr = foptr;
            indexx = 0; // the selected polygon is the first one
        }
        break;
    case 27: // <ESC>
        exit(0);
        break;
    default:
        printf("%d %c\n", key, key);
    }

    // The screen must be drawn to show the new triangle
    glutPostRedisplay();
}

int main(int argc, char **argv)
{
    int retval;

    printf(" Triangeluak: barneko puntuak eta testura\n Triángulos con puntos internos y textura \n");
    printf("Press <ESC> to finish\n");
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("KBG/GO praktika");

    glutDisplayFunc(marraztu);
    glutKeyboardFunc(teklatua);
    /* we put the information of the texture in the buffer pointed by bufferra. The dimensions of the texture are loaded into dimx and dimy */
    retval = load_ppm("testura.ppm", &bufferra, &dimx, &dimy);
    if (retval != 1)
    {
        printf("Ez dago texturaren fitxategia (testura.ppm)\n");
        exit(-1);
    }

    glClearColor(0.0f, 0.0f, 0.7f, 1.0f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glEnable(GL_DEPTH_TEST); // activar el test de profundidad (Z-buffer)
    denak = 1;
    lineak = 1;
    objektuak = 1;
    foptr = 0;
    sel_ptr = 0;
    aldaketa = 'r';
    ald_lokala = 1;

    // TODO: Temporal, cargar una camara como un objeto. Buscar otra manera mas simple
    read_from_file("camara.txt");
    camara_ptr = sel_ptr;

    if (argc > 1)
        read_from_file(argv[1]);
    else
        read_from_file("adibideak.txt");

    glutMainLoop();

    return 0;
}
