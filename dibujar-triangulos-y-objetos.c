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
// #include "cargar-triangulo.h"
#include "obj.h"

#define DESPLAZAMIENTO_TRANSLACION 10
#define ANGULO_ROTACION 3.14159265358979323846 / 32
#define PROPORCION_ESCALADO 1.2
#define DISTANCIA_MINIMA_ANALISIS 30

#define CAMARA_CONFIG_NEAR 5.0
#define CAMARA_CONFIG_FAR 1000.0
#define CAMARA_CONFIG_LEFT -5.0
#define CAMARA_CONFIG_RIGHT 5.0
#define CAMARA_CONFIG_TOP 5.0
#define CAMARA_CONFIG_BOTTOM -5.0

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

#define LISTA_OBJETOS 0
#define LISTA_CAMARAS 1
// #define LISTA_LUCES 2

#define CAMARA_PERSPECTIVA 0
#define CAMARA_PARALELA 1

#define CAMARA_MOD_VUELO 0
#define CAMARA_MOD_ANALISIS 1

// typedef struct mlist
// {
//     double m[16];
//     struct mlist *hptr;
// } mlist;

/*
typedef struct triobj
{
    hiruki *triptr;
    int num_triangles;
    mlist *mptr;
    struct triobj *hptr;
    unsigned char *color;
} triobj;
*/

// testuraren informazioa
// información de textura

extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int *dimyptr);
unsigned char *bufferra;
int dimx, dimy;

int indexx;
// hiruki *triangulosptr;
mlist *mcsr_ptr;
mlist *mmodelview_ptr;
mlist *mperspectiva_ptr;

object3d *objetosptr; // Lista de objetos ( Apunta al primer objeto )
object3d *obj_ptr;    // Puntero al objeto seleccionado de la lista

object3d *camarasptr; // Lista de camaras ( Apunta a la primera camara )
object3d *camara_ptr; // Puntero a la camara seleccionada de la lista

int lista_activa; // Indicador de la lista actualmente activa ( Camaras, objetos, ... ). Solo lectura

object3d **foptr;   // Puntero al puntero del primer objeto de la lista activa (independiente de la lista)
object3d **sel_ptr; // Puntero al puntero del objeto actualmente seleccionado (independiente de la lista)

int camara_activa; // Si estamos viendo desde la camara o no ( en caso negativo vemos desde objeto )
int tipo_camara;
int modo_camara;
int denak;
int lineak;
int objektuak;
char aldaketa;
int ald_lokala;

char fitxiz[100];

int ultimo_es_visible;
int dibujar_no_visible;
int dibujar_normales;

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

void dibujar_linea(punto p1, punto p2, color3 color)
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

    if (pcortemayor->x - pcortemenor->x > 1000)
        return;
    if (abs(pcortemayor->y - pcortemenor->y) > 1000)
        return;
    if (abs(pcortemayor->z - pcortemenor->z) > 1000)
        return;

    for (j = 1; j > 0; j -= cambioj)
    {
        pcalculado.x = j * pcortemayor->x + (1 - j) * pcortemenor->x;
        pcalculado.y = j * pcortemayor->y + (1 - j) * pcortemenor->y;
        pcalculado.z = j * pcortemayor->z + (1 - j) * pcortemenor->z;
        pcalculado.u = j * pcortemayor->u + (1 - j) * pcortemenor->u;
        pcalculado.v = j * pcortemayor->v + (1 - j) * pcortemenor->v;

        // TODO: Por ahora me vale para mejorar un poco el rendimiento
        if (abs(pcalculado.x) > 500 || abs(pcalculado.y) > 500)
            continue;

        glBegin(GL_POINTS);
        if (ultimo_es_visible == 1)
        {
            r = color.r;
            g = color.g;
            b = color.b;
        }
        else
        {
            r = 255;
            g = 0;
            b = 0;
        }
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
        printf("%lf, %lf, %lf, %lf\n", (*sel_ptr)->mptr->m[i * 4], (*sel_ptr)->mptr->m[i * 4 + 1], (*sel_ptr)->mptr->m[i * 4 + 2],
               (*sel_ptr)->mptr->m[i * 4 + 3]);
}

void print_matrizea2(char *str, mlist *matrizea)
{
    int i;

    printf("%s\n", str);
    for (i = 0; i < 4; i++)
        printf("%lf, %lf, %lf, %lf\n", matrizea->m[i * 4], matrizea->m[i * 4 + 1], matrizea->m[i * 4 + 2],
               matrizea->m[i * 4 + 3]);
}

void print_estado()
{
    printf("\n");
    printf("------------(COMIENZO ESTADO)------------\n");

    print_matrizea2("Posición cámara:", camara_ptr->mptr);
    printf("\n");

    printf("Modo cámara (g): ");
    if (modo_camara == CAMARA_MOD_VUELO)
        printf("VUELO\n");
    else
        printf("ANALISIS\n");

    printf("Tipo camara (p): ");
    if (tipo_camara == CAMARA_PERSPECTIVA)
        printf("PERSPECTIVA\n");
    else
        printf("PARALELO\n");

    printf("Lista activa (c): ");
    if (lista_activa == LISTA_CAMARAS)
        printf("CAMARAS\n");
    else
        printf("OBJETOS\n");

    printf("Dibujando caras traseras (b): ");
    if (dibujar_no_visible == 1)
        printf("SI\n");
    else
        printf("NO\n");

    printf("Dibujando normales (n): ");
    if (dibujar_normales == 1)
        printf("SI\n");
    else
        printf("NO\n");

    printf("------------(FIN ESTADO)------------\n");
}

void mxp(punto *pptr, double m[16], punto p)
{
    pptr->x = m[0] * p.x + m[1] * p.y + m[2] * p.z + m[3];
    pptr->y = m[4] * p.x + m[5] * p.y + m[6] * p.z + m[7];
    pptr->z = m[8] * p.x + m[9] * p.y + m[10] * p.z + m[11];
    pptr->u = p.u;
    pptr->v = p.v;
}

void mxv(punto *pptr, double m[16], vertex v)
{
    pptr->x = m[0] * v.coord.x + m[1] * v.coord.y + m[2] * v.coord.z + m[3];
    pptr->y = m[4] * v.coord.x + m[5] * v.coord.y + m[6] * v.coord.z + m[7];
    pptr->z = m[8] * v.coord.x + m[9] * v.coord.y + m[10] * v.coord.z + m[11];
    pptr->u = v.u;
    pptr->v = v.v;
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

void calcular_mcsr(mlist *mresptr, double mCamara[16])
{
    int i;
    for (i = 0; i < 16; i++)
        mresptr->m[i] = 0;

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
    mresptr->m[0] = mCamara[0];
    mresptr->m[1] = mCamara[4];
    mresptr->m[2] = mCamara[8];

    // Columna Yc
    mresptr->m[4] = mCamara[1];
    mresptr->m[5] = mCamara[5];
    mresptr->m[6] = mCamara[9];

    // Columna Zc
    mresptr->m[8] = mCamara[2];
    mresptr->m[9] = mCamara[6];
    mresptr->m[10] = mCamara[10];

    // Columna -(X,Y,Z)c * Pc
    mresptr->m[3] = -(mCamara[0] * mCamara[3] + mCamara[4] * mCamara[7] + mCamara[8] * mCamara[11]);   // -Xc * Pc
    mresptr->m[7] = -(mCamara[1] * mCamara[3] + mCamara[5] * mCamara[7] + mCamara[9] * mCamara[11]);   // -Yc * Pc
    mresptr->m[11] = -(mCamara[2] * mCamara[3] + mCamara[6] * mCamara[7] + mCamara[10] * mCamara[11]); // -Zc * Pc

    mresptr->m[15] = 1;
}

void calcular_mmodelview(mlist *mresptr, double mcsr[16], double mobj[16])
{
    int i;
    for (i = 0; i < 16; i++)
        mresptr->m[i] = 0;
    mxm(mresptr->m, mcsr, mobj);
}

void calcular_mperspectiva(mlist *mresptr, float n, float f, float r, float l, float t, float b)
{
    int i;
    float r_menos_l, t_menos_b, f_menos_n;

    for (i = 0; i < 16; i++)
        mresptr->m[i] = 0;

    r_menos_l = r - l;
    t_menos_b = t - b;
    f_menos_n = f - n;

    if (r_menos_l == 0)
    {
        printf("r_menos_l == 0\n");
        return;
    }
    if (t_menos_b == 0)
    {
        printf("t_menos_b == 0\n");
        return;
    }
    if (f_menos_n == 0)
    {
        printf("f_menos_n == 0\n");
        return;
    }

    mresptr->m[0] = (2 * n) / r_menos_l;
    mresptr->m[2] = (r + l) / r_menos_l;
    mresptr->m[5] = (2 * n) / t_menos_b;
    mresptr->m[6] = (t + b) / t_menos_b;
    mresptr->m[10] = -(f + n) / f_menos_n;
    mresptr->m[11] = (-2 * f * n) / f_menos_n;
    mresptr->m[14] = -1;

    // print_matrizea2("Mpers\n",mresptr);
}

// Devuelve 0 si el punto está dentro del campo de vision
// Devuelve -1 si la cuarta componente es 0 y la division da infinito
// Devuelve -2 si el punto esta detras de la camara
int aplicar_mperspectiva(punto *pptr, double m[16])
{

    punto ptemp;
    float w;

    /*
        m[16] posiciones
        [ 0  1  2  3  |
        | 4  5  6  7  |
        | 8  9  10 11 |
        | 12 13 14 15 ]
    */

    ptemp.x = m[0] * pptr->x + m[1] * pptr->y + m[2] * pptr->z + m[3];
    ptemp.y = m[4] * pptr->x + m[5] * pptr->y + m[6] * pptr->z + m[7];
    ptemp.z = m[8] * pptr->x + m[9] * pptr->y + m[10] * pptr->z + m[11];
    // w = m[12] * pptr->x + m[13] * pptr->y + m[14] * pptr->z + m[15];
    w = m[14] * pptr->z;
    ptemp.u = pptr->u;
    ptemp.v = pptr->v;

    if (abs(w) < 1)
    {
        // printf("INTERSECCION CON PLANO NEAR\n");
        return -1; // Interseccion con plano near
    }

    // printf("x=%f z=%f w=%f ptemp (%.3f,%.3f,%.3f)\n", pptr->x, pptr->z, w, ptemp.x / w, ptemp.y / w, -ptemp.z / w);

    // printf(" 1 ptemp.x %f\n", ptemp.x);
    // printf(" 1 pptr->x %f\n", pptr->x);
    // printf(" 1 w %f\n", w);
    pptr->x = (ptemp.x / w) * 500;
    pptr->y = (ptemp.y / w) * 500;
    pptr->z = -(ptemp.z / abs(w));
    pptr->u = ptemp.u;
    pptr->v = ptemp.v;

    // z visible [-500,-CAMARA_CONFIG_NEAR)
    if (pptr->z >= 0 || pptr->z <= -1)
    {
        return -2;
    }

    pptr->z *= 500;

    return 0;

    // printf(" 2 pptr->x %f\n", pptr->z);
}

// TODO: si miras a un objeto en tu misma pos va a dar error, siendo los resultados nan
// Arreglo temporal: Poner la identidad en los vectores
void look_at(object3d *observadorptr, object3d *objetivoptr)
{
    int i;
    mlist *nueva_matrizptr = (mlist *)malloc(sizeof(mlist));

    for (i = 0; i < 16; i++)
        nueva_matrizptr->m[i] = 0;

    // Parametros: E, At, (Vp = (0,1,0))
    double e[3] = {observadorptr->mptr->m[3], observadorptr->mptr->m[7], observadorptr->mptr->m[11]};
    double at[3] = {objetivoptr->mptr->m[3], objetivoptr->mptr->m[7], objetivoptr->mptr->m[11]};
    double vp[3] = {0, 1, 0};
    // Buscamos: Zc, Xc, Yc
    double z_c[3];
    double x_c[3];
    double y_c[3];

    // Zc = (E - At) / || E - At ||
    double e_at[3] = {e[0] - at[0], e[1] - at[1], e[2] - at[2]};
    double at_e[3] = {at[0] - e[0], at[1] - e[1], at[2] - e[2]};
    double mod_e_at = sqrt(pow(e_at[0], 2) + pow(e_at[1], 2) + pow(e_at[2], 2));

    // At - E
    // E - At /= Vp
    if (at_e[0] == vp[0] && at_e[1] == vp[1] && at_e[2] == vp[2])
        printf("IGUALES   AT-E=VP\n");
    if (e_at[0] == vp[0] && e_at[1] == vp[1] && e_at[2] == vp[2])
        printf("IGUALES   E-AT=VP\n");

    if (mod_e_at == 0)
    {
        z_c[0] = 0;
        z_c[1] = 0;
        z_c[2] = 1;
    }
    else
    {
        z_c[0] = e_at[0] / mod_e_at;
        z_c[1] = e_at[1] / mod_e_at;
        z_c[2] = e_at[2] / mod_e_at;
    }

    // Xc = (Vp /\ Zc) / ||(Vp /\ Zc)||
    double vp_zc[3] = {vp[1] * z_c[2] - vp[2] * z_c[1], -(vp[0] * z_c[2] - vp[2] * z_c[0]), vp[0] * z_c[1] - vp[1] * z_c[0]};
    double mod_vp_zc = sqrt(pow(vp_zc[0], 2) + pow(vp_zc[1], 2) + pow(vp_zc[2], 2));
    if (mod_vp_zc == 0)
    {
        x_c[0] = 1;
        x_c[1] = 0;
        x_c[2] = 0;
    }
    else
    {
        x_c[0] = vp_zc[0] / mod_vp_zc;
        x_c[1] = vp_zc[1] / mod_vp_zc;
        x_c[2] = vp_zc[2] / mod_vp_zc;
    }

    // Yc = Zc /\ Xc
    double zc_xc[3] = {z_c[1] * x_c[2] - z_c[2] * x_c[1], -(z_c[0] * x_c[2] - z_c[2] * x_c[0]), z_c[0] * x_c[1] - z_c[1] * x_c[0]};
    y_c[0] = zc_xc[0];
    y_c[1] = zc_xc[1];
    y_c[2] = zc_xc[2];

    // Ahora que tenemos Xc, Yc y Zc, construimos la nueva matriz de la camara
    for (i = 0; i < 3; i++)
        nueva_matrizptr->m[i * 4] = x_c[i];
    for (i = 0; i < 3; i++)
        nueva_matrizptr->m[(i * 4) + 1] = y_c[i];
    for (i = 0; i < 3; i++)
        nueva_matrizptr->m[(i * 4) + 2] = z_c[i];
    for (i = 0; i < 3; i++)
        nueva_matrizptr->m[(i * 4) + 3] = e[i];

    nueva_matrizptr->m[15] = 1;

    nueva_matrizptr->hptr = observadorptr->mptr;
    observadorptr->mptr = nueva_matrizptr;
}

// Teniendo el cuenta la lista activa cambia sel_ptr al siguiente elemento
void siguiente_elemento_lista()
{
    if ((*foptr) != 0) // objekturik gabe ez du ezer egin behar
                       // si no hay objeto no hace nada
    {
        (*sel_ptr) = (*sel_ptr)->hptr;
        /*The selection is circular, thus if we move out of the list we go back to the first element*/
        if ((*sel_ptr) == 0)
            (*sel_ptr) = (*foptr);
        indexx = 0; // the selected polygon is the first one
    }
}

// Cambia cual es la lista activa, haciendo que sel_ptr apunte al elemento seleccionado de esa lista
// Y foptr al primer elemento de esa lista
void cambiar_lista_activa(int lista)
{
    char *nombre_lista;

    switch (lista)
    {
    case LISTA_OBJETOS:
        sel_ptr = &obj_ptr;
        foptr = &objetosptr;
        nombre_lista = "LISTA_OBJETOS";
        break;
    case LISTA_CAMARAS:
        sel_ptr = &camara_ptr;
        foptr = &camarasptr;
        nombre_lista = "LISTA_CAMARAS";
        break;
    }

    lista_activa = lista;

    printf("Se ha cambiado la lista activa a: %s\n", nombre_lista);
}

int es_visible(object3d *optr, int i)
{
    face *fptr;
    object3d *observadorptr;
    mlist matriz_csr_objeto;
    mlist matriz_observador;

    if (camara_activa == 1)
        observadorptr = camara_ptr;
    else
        observadorptr = obj_ptr;

    if (i >= optr->num_faces)
        return 0;
    fptr = optr->face_table + i;

    if (fptr->num_vertices == 0 || optr->num_vertices == 0)
        return 0;

    calcular_mcsr(&matriz_csr_objeto, optr->mptr->m);
    mxm(matriz_observador.m, matriz_csr_objeto.m, observadorptr->mptr->m);

    // if (v * n > 0) dibujar
    // else no dibujar

    if (tipo_camara == CAMARA_PERSPECTIVA)
    {
        double v[3] = {matriz_observador.m[3] - optr->vertex_table[fptr->vertex_ind_table[0]].coord.x, matriz_observador.m[7] - optr->vertex_table[fptr->vertex_ind_table[0]].coord.y, matriz_observador.m[11] - optr->vertex_table[fptr->vertex_ind_table[0]].coord.z}; // observador - punto (0,0,0)
        double v_n = v[0] * fptr->N[0] + v[1] * fptr->N[1] + v[2] * fptr->N[2];                                                                                                                                                                                          // v * n
        return v_n > 0;
    }
    else // CAMARA_PARALELA
    {
        // En paralelo unicamente nos interesa hacia donde esta mirando la camara, es decir, su vector de direccion hacia delante ( la z )
        double v_n = matriz_observador.m[2] * fptr->N[0] + matriz_observador.m[6] * fptr->N[1] + matriz_observador.m[10] * fptr->N[2]; // vector_z_camara * n
        return v_n > 0;
    }
}

void dibujar_triangulo(object3d *optr, int i)
{
    face *fptr;

    punto *pgoiptr, *pbeheptr, *perdiptr;
    punto p1, p2, p3;
    punto p2_vnormal;
    punto p2_vnormal_transformado;

    float t, s;
    float cambiot, cambios;

    punto pcorte1, pcorte2;

    int res_mpers_vnormal;

    if (i >= optr->num_faces)
        return;
    fptr = optr->face_table + i;

    if (fptr->num_vertices == 0 || optr->num_vertices == 0)
        return;

    ultimo_es_visible = es_visible(optr, i);

    if (ultimo_es_visible != 1 && dibujar_no_visible != 1)
        return;

    // (Mcsr * Mobj) * Obj
    mxv(&p1, mmodelview_ptr->m, optr->vertex_table[fptr->vertex_ind_table[0]]);
    mxv(&p2, mmodelview_ptr->m, optr->vertex_table[fptr->vertex_ind_table[1]]);
    mxv(&p3, mmodelview_ptr->m, optr->vertex_table[fptr->vertex_ind_table[2]]);

    // Si la camara esta en perspectiva, aplicar (Mp * Mcsr * Mobj * Obj)
    if (tipo_camara == CAMARA_PERSPECTIVA)
    {

        if (aplicar_mperspectiva(&p1, mperspectiva_ptr->m) != 0)
            return;
        if (aplicar_mperspectiva(&p2, mperspectiva_ptr->m) != 0)
            return;
        if (aplicar_mperspectiva(&p3, mperspectiva_ptr->m) != 0)
            return;
    }

    if (lineak == 1)
    {
        glBegin(GL_POLYGON);

        if (ultimo_es_visible != 1)
            glColor3ub(255, 0, 0);
        else
            glColor3ub(255, 255, 255);

        glVertex3d(p1.x, p1.y, p1.z);
        glVertex3d(p2.x, p2.y, p2.z);
        glVertex3d(p3.x, p3.y, p3.z);

        // printf("p1 (%f,%f,%f) p2 (%f,%f,%f) p3 (%f,%f,%f)\n",p1.x,p1.y,p1.z,p2.x,p2.y,p2.z,p3.x,p3.y,p3.z);

        glEnd();

        if (dibujar_normales == 1)
        {

            // NORMALES DE LAS CARAS
            p2_vnormal.x = optr->vertex_table[fptr->vertex_ind_table[0]].coord.x + fptr->N[0] * 20;
            p2_vnormal.y = optr->vertex_table[fptr->vertex_ind_table[0]].coord.y + fptr->N[1] * 20;
            p2_vnormal.z = optr->vertex_table[fptr->vertex_ind_table[0]].coord.z + fptr->N[2] * 20;
            p2_vnormal.u = 0;
            p2_vnormal.v = 0;
            mxp(&p2_vnormal_transformado, mmodelview_ptr->m, p2_vnormal);
            if (tipo_camara == CAMARA_PERSPECTIVA)
            {
                res_mpers_vnormal = aplicar_mperspectiva(&p2_vnormal_transformado, mperspectiva_ptr->m);
                // if (p2_vnormal_transformado.z < -1 || p2_vnormal_transformado.z >= 0)
                //     res_mpers_vnormal = -1;
            }
            else
            {
                res_mpers_vnormal = 0;
            }

            if (res_mpers_vnormal == 0)
            {
                glBegin(GL_LINES);
                glVertex3d(p1.x, p1.y, p1.z);
                glVertex3d(p2_vnormal_transformado.x, p2_vnormal_transformado.y, p2_vnormal_transformado.z);
                glEnd();
            }

            // NORMALES DE LOS VERTICES
            for (i = 0; i < fptr->num_vertices; i++)
            {
                // printf("optr->vertex_table[fptr->vertex_ind_table[i]].coord.x=%f\n",optr->vertex_table[fptr->vertex_ind_table[i]].coord.x);
                // printf("optr->vertex_table[fptr->vertex_ind_table[i]].N[0]=%f\n",optr->vertex_table[fptr->vertex_ind_table[i]].N[0]);
                p2_vnormal.x = optr->vertex_table[fptr->vertex_ind_table[i]].coord.x + optr->vertex_table[fptr->vertex_ind_table[i]].N[0] * 3;
                p2_vnormal.y = optr->vertex_table[fptr->vertex_ind_table[i]].coord.y + optr->vertex_table[fptr->vertex_ind_table[i]].N[1] * 3;
                p2_vnormal.z = optr->vertex_table[fptr->vertex_ind_table[i]].coord.z + optr->vertex_table[fptr->vertex_ind_table[i]].N[2] * 3;
                p2_vnormal.u = 0;
                p2_vnormal.v = 0;

                mxp(&p2_vnormal_transformado, mmodelview_ptr->m, p2_vnormal);
                if (tipo_camara == CAMARA_PERSPECTIVA)
                {
                    res_mpers_vnormal = aplicar_mperspectiva(&p2_vnormal_transformado, mperspectiva_ptr->m);
                    // if (p2_vnormal_transformado.z < -1 || p2_vnormal_transformado.z >= 0)
                    //     res_mpers_vnormal = -1;
                }
                else
                {
                    res_mpers_vnormal = 0;
                }

                if (res_mpers_vnormal == 0)
                {
                    glBegin(GL_LINES);
                    glColor3ub(134,134,134);
                    glVertex3d(p1.x, p1.y, p1.z);
                    glVertex3d(p2_vnormal_transformado.x, p2_vnormal_transformado.y, p2_vnormal_transformado.z);
                    glEnd();
                }
            }
        }

        return;
    }
    //  else

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

        dibujar_linea(pcorte1, pcorte2, optr->rgb);
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

        dibujar_linea(pcorte1, pcorte2, optr->rgb);
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

        dibujar_linea(pcorte1, pcorte2, optr->rgb);
    }
}

static void marraztu(void)
{
    float u, v;
    int i, j;
    mlist mcsr;
    mlist mmodelview;
    mlist mperspectiva;
    mlist mtemp;
    object3d *auxptr;
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
    // glOrtho(-500.0, 500.0, -500.0, 500.0, -500.0, 500.0);
    // glOrtho(-500.0, 500.0, -500.0, 500.0, 0.0, 500.0); // Para no dibujar lo que está detrás
    // glOrtho(CAMARA_CONFIG_LEFT, CAMARA_CONFIG_RIGHT, CAMARA_CONFIG_BOTTOM, CAMARA_CONFIG_TOP, CAMARA_CONFIG_NEAR, CAMARA_CONFIG_FAR); // Para no dibujar lo que está detrás
    if (tipo_camara == CAMARA_PARALELA)
    {
        glOrtho(-500.0, 500.0, -500.0, 500.0, 0, 500.0);
    }
    else
        glOrtho(-500.0, 500.0, -500.0, 500.0, -500.0, 500.0);

    mcsr_ptr = &mcsr;
    mmodelview_ptr = &mmodelview;

    // Si tenemos camara multiplicar la 'Mcsr' a la matriz del objeto 'Mobj'
    if (camara_ptr != 0 && camara_activa == 1)
    {
        calcular_mcsr(mcsr_ptr, camara_ptr->mptr->m);
    }
    else
    {
        calcular_mcsr(mcsr_ptr, obj_ptr->mptr->m);
    }

    // triangulosptr = (*sel_ptr)->triptr;
    if (objektuak == 1)
    {
        if (denak == 1)
        {
            for (auxptr = objetosptr; auxptr != 0; auxptr = auxptr->hptr)
            {
                calcular_mmodelview(mmodelview_ptr, mcsr_ptr->m, auxptr->mptr->m);
                for (i = 0; i < auxptr->num_faces; i++)
                {
                    dibujar_triangulo(auxptr, i);
                }
            }

            for (auxptr = camarasptr; auxptr != 0; auxptr = auxptr->hptr)
            {
                calcular_mmodelview(mmodelview_ptr, mcsr_ptr->m, auxptr->mptr->m);
                for (i = 0; i < auxptr->num_faces; i++)
                {
                    dibujar_triangulo(auxptr, i);
                }

                // print_matrizea2("Mcamara\n", camara_ptr->mptr);
            }
        }
        else
        {
            calcular_mmodelview(mmodelview_ptr, mcsr_ptr->m, (*sel_ptr)->mptr->m);
            for (i = 0; i < (*sel_ptr)->num_faces; i++)
            {
                dibujar_triangulo((*sel_ptr), i);
            }
        }
    }
    else
    {
        dibujar_triangulo((*sel_ptr), indexx);
    }
    glFlush();
}

void read_from_file(char *fitx, int tipo_lista)
{
    int i, retval;
    object3d *optr;

    cambiar_lista_activa(tipo_lista);
    // Se va a añadir el objeto nuevo a la lista que toca
    // Ademas se va a cambiar el objeto seleccionado al nuevo
    // haciendo que las transformaciones se apliquen sobre este.
    // Esto no hace que si estabas viendo desde la camara cambie
    // a ver desde el objeto o al reves.

    // printf("%s fitxategitik datuak hartzera\n",fitx);
    optr = (object3d *)malloc(sizeof(object3d));
    // retval = cargar_triangulos_color(fitx, &(optr->num_faces), &(optr->triptr), &(optr->color));
    retval = read_wavefront(fitx, optr);
    if (retval != 0)
    {
        printf("%s fitxategitik datuak hartzerakoan arazoak izan ditut\n    Problemas al leer\n cod:%d\n", fitxiz, retval);
        free(optr);
        exit(-1);
    }
    else
    {
        // triangulosptr = optr->triptr;
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
        // printf("objektu zerrendara doa informazioa...\n");ultimo_listaptrptr
        optr->hptr = (*foptr); // el siguiente al objeto es el que antes era el primero
        (*foptr) = optr;       // foptr(cambiado a un ptrptr que apunta dependiendo de la lista) apunta al primer objeto ( el ultimo cargado )
        (*sel_ptr) = optr;
        if (retval == 9)
            printf("COLOR (%d,%d,%d)\n", optr->rgb.r, optr->rgb.g, optr->rgb.b);
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

void rotacion_respecto_punto(mlist *matriz_transformacion, object3d *pA, int eje, int dir, double angulo)
{
    int i;
    float angulo_x_dir;
    mlist m_translacion;
    mlist m_deshacer_translacion;
    mlist m_rotacion_eje;
    mlist m_temp;
    double eje_rotacion[3];
    double cos_angulo, sin_angulo;

    for (i = 0; i < 16; i++)
    {
        m_temp.m[i] = 0;
        m_rotacion_eje.m[i] = 0;
        matriz_transformacion->m[i] = 0;
    }

    if (dir == DIR_ADELANTE)
        angulo_x_dir = angulo;
    else // dir == DIR_ATRAS
        angulo_x_dir = (2 * 3.14159) - angulo;

    if (eje == EJE_X)
    {
        eje_rotacion[0] = (*sel_ptr)->mptr->m[0];
        eje_rotacion[1] = (*sel_ptr)->mptr->m[4];
        eje_rotacion[2] = (*sel_ptr)->mptr->m[8];
    }
    else if (eje == EJE_Y)
    {
        eje_rotacion[0] = (*sel_ptr)->mptr->m[1];
        eje_rotacion[1] = (*sel_ptr)->mptr->m[5];
        eje_rotacion[2] = (*sel_ptr)->mptr->m[9];
    }
    else // EJE_Z
    {
        eje_rotacion[0] = (*sel_ptr)->mptr->m[2];
        eje_rotacion[1] = (*sel_ptr)->mptr->m[6];
        eje_rotacion[2] = (*sel_ptr)->mptr->m[10];
    }

    translacion(&m_translacion, EJE_X, DIR_ADELANTE, 0);
    translacion(&m_deshacer_translacion, EJE_X, DIR_ADELANTE, 0);

    // Hacer que el objeto sea el origen
    m_translacion.m[3] = -pA->mptr->m[3];   // x
    m_translacion.m[7] = -pA->mptr->m[7];   // y
    m_translacion.m[11] = -pA->mptr->m[11]; // z

    // Deshacer la translacion
    m_deshacer_translacion.m[3] = pA->mptr->m[3];   // x
    m_deshacer_translacion.m[7] = pA->mptr->m[7];   // y
    m_deshacer_translacion.m[11] = pA->mptr->m[11]; // z

    // Rotar respecto a un punto ( siendo el punto el punto de atencion )
    cos_angulo = cos(angulo_x_dir);
    sin_angulo = sin(angulo_x_dir);

    m_rotacion_eje.m[0] = cos_angulo + (1 - cos_angulo) * pow(eje_rotacion[0], 2);
    m_rotacion_eje.m[1] = (1 - cos_angulo) * eje_rotacion[0] * eje_rotacion[1] - eje_rotacion[2] * sin_angulo;
    m_rotacion_eje.m[2] = (1 - cos_angulo) * eje_rotacion[0] * eje_rotacion[2] + eje_rotacion[1] * sin_angulo;
    m_rotacion_eje.m[4] = (1 - cos_angulo) * eje_rotacion[0] * eje_rotacion[1] + eje_rotacion[2] * sin_angulo;
    m_rotacion_eje.m[5] = cos_angulo + (1 - cos_angulo) * pow(eje_rotacion[1], 2);
    m_rotacion_eje.m[6] = (1 - cos_angulo) * eje_rotacion[1] * eje_rotacion[2] - eje_rotacion[0] * sin_angulo;
    m_rotacion_eje.m[8] = (1 - cos_angulo) * eje_rotacion[0] * eje_rotacion[2] - eje_rotacion[1] * sin_angulo;
    m_rotacion_eje.m[9] = (1 - cos_angulo) * eje_rotacion[1] * eje_rotacion[2] + eje_rotacion[0] * sin_angulo;
    m_rotacion_eje.m[10] = cos_angulo + (1 - cos_angulo) * pow(eje_rotacion[2], 2);
    m_rotacion_eje.m[15] = 1;

    mxm(m_temp.m, m_rotacion_eje.m, m_translacion.m);
    mxm(matriz_transformacion->m, m_deshacer_translacion.m, m_temp.m);

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
        mxm(nueva_matrizptr->m, (*sel_ptr)->mptr->m, matriz_transformacionptr->m);
    }
    else // SISTEMA_MUNDO
    {
        // Multiplicar por la izquierda
        mxm(nueva_matrizptr->m, matriz_transformacionptr->m, (*sel_ptr)->mptr->m);
    }

    nueva_matrizptr->hptr = (*sel_ptr)->mptr;
    (*sel_ptr)->mptr = nueva_matrizptr;
}

// Hace que la camara no pueda estar a menos de X distancia del objetivo
// ajustando la z para que sea asi
void ajustar_distancia_analisis()
{
    mlist matriz_transformacion;
    double distancia;
    double pos_cam[3] = {camara_ptr->mptr->m[3], camara_ptr->mptr->m[7], camara_ptr->mptr->m[11]};
    double pos_obj[3] = {obj_ptr->mptr->m[3], obj_ptr->mptr->m[7], obj_ptr->mptr->m[11]};

    distancia = sqrt(pow(pos_obj[0] - pos_cam[0], 2) + pow(pos_obj[1] - pos_cam[1], 2) + pow(pos_obj[2] - pos_cam[2], 2));

    if (distancia < DISTANCIA_MINIMA_ANALISIS)
    {

        translacion(&matriz_transformacion, EJE_Z, DIR_ADELANTE, (DISTANCIA_MINIMA_ANALISIS - distancia));
        aplicar_transformacion(&matriz_transformacion, SISTEMA_LOCAL);
    }

    // printf("DISTANCIA:%f\n", distancia);
}

void tratar_transformacion_modo_analisis(int eje, int dir)
{
    mlist matriz_transformacion;
    int sistema_ref;
    // Si estamos con camara_activa y en modo analisis
    // solo disponemos de dos transformaciones ( y limitadas ):
    // 1 . Rotacion :  Eje X y Eje Y
    // 2 . Translacion : Eje Z

    switch (aldaketa)
    {
    case TRANSLACION:
        if (eje != EJE_Z)
            return;

        translacion(&matriz_transformacion, eje, dir, DESPLAZAMIENTO_TRANSLACION);
        sistema_ref = SISTEMA_LOCAL;
        break;
    case ROTACION:
        if (eje == EJE_X)
            rotacion_respecto_punto(&matriz_transformacion, obj_ptr, EJE_Y, dir, ANGULO_ROTACION);
        else if (eje == EJE_Y)
            rotacion_respecto_punto(&matriz_transformacion, obj_ptr, EJE_X, dir, ANGULO_ROTACION);
        else
            return;
        sistema_ref = SISTEMA_MUNDO;
        break;
    default:
        return;
    }

    aplicar_transformacion(&matriz_transformacion, sistema_ref);
    // print_matrizea2("Mcam:\n", camara_ptr->mptr);
    // print_matrizea2("Mobj:\n", obj_ptr->mptr);
    ajustar_distancia_analisis();
    print_estado();
}

void tratar_transformacion(int eje, int dir)
{
    mlist matriz_transformacion;

    // Si estamos con camara_activa y en modo analisis,
    // vamos a separar la logica para las transformaciones en este caso

    if (camara_activa == 1 && modo_camara == CAMARA_MOD_ANALISIS && lista_activa == LISTA_CAMARAS)
    {
        tratar_transformacion_modo_analisis(eje, dir);
        return;
    }

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
        return;
    }

    // Las transformaciones a las camaras siempre se hacen en el sistema local,
    // independientemente del sistema de referencia activo
    if (lista_activa == LISTA_CAMARAS)
        aplicar_transformacion(&matriz_transformacion, SISTEMA_LOCAL); // La camara siempre en sis.local
    else
        aplicar_transformacion(&matriz_transformacion, ald_lokala);

    print_estado();
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

    if ((*sel_ptr)->mptr->hptr != 0)
    {
        matriz_a_borrarptr = (*sel_ptr)->mptr;
        (*sel_ptr)->mptr = (*sel_ptr)->mptr->hptr;
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
        if ((*foptr) != 0) // objekturik ez badago ezer ez du egin behar
                           // si no hay objeto que no haga nada
        {
            indexx++; // azkena bada lehenengoa bihurtu
                      // pero si es el último? hay que controlarlo!
            if (indexx == (*sel_ptr)->num_faces)
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
    case 'G':
    case 'g':
        if (lista_activa == LISTA_CAMARAS) // ANALISIS <-> VUELO
        {
            if (modo_camara == CAMARA_MOD_VUELO)
                modo_camara = CAMARA_MOD_ANALISIS;
            else
                modo_camara = CAMARA_MOD_VUELO;

            if (modo_camara == CAMARA_MOD_VUELO)
                printf("Cambiado a modo vuelo\n");
            else
            {
                look_at(camara_ptr, obj_ptr);
                ajustar_distancia_analisis();
                printf("Cambiado a modo analisis\n");
            }
        }
        else
        { // S.REF.LOCAL <-> S.REF.MUNDO

            if (ald_lokala == 1)
                ald_lokala = 0;
            else
                ald_lokala = 1;

            if (ald_lokala == SISTEMA_LOCAL)
                printf("Cambiado a sistema local\n");
            else
                printf("Cambiado al sistema del mundo\n");
        }

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
    case 'c':
        if (lista_activa == LISTA_CAMARAS)
        {
            cambiar_lista_activa(LISTA_OBJETOS);
            if (modo_camara == CAMARA_MOD_ANALISIS)
            {
                modo_camara = CAMARA_MOD_VUELO;
            }
        }
        else if (camara_activa == 1)
            cambiar_lista_activa(LISTA_CAMARAS);
        break;
    case 'C':
        if (camara_activa == 0)
        {
            printf("Viendo el mundo desde la cámara\n");
            camara_activa = 1;
        }
        else
        {
            printf("Viendo el mundo desde el objeto seleccionado\n");
            camara_activa = 0;
            cambiar_lista_activa(LISTA_OBJETOS);
        }
        break;
    case 'P':
    case 'p':
        if (tipo_camara == CAMARA_PARALELA)
        {
            printf("Modo camara: CAMARA_PERSPECTIVA\n");
            tipo_camara = CAMARA_PERSPECTIVA;
        }
        else
        {
            printf("Modo camara: CAMARA_PARALELA\n");
            tipo_camara = CAMARA_PARALELA;
        }
        break;
    case 'b':
        if (dibujar_no_visible == 1)
            dibujar_no_visible = 0;
        else
            dibujar_no_visible = 1;
        break;
    case 'n':
        if (dibujar_normales == 1)
            dibujar_normales = 0;
        else
            dibujar_normales = 1;
        break;
    case 'f':
        /*Ask for file*/
        printf("idatzi fitxategi izena\n");
        scanf("%s", &(fitxiz[0]));
        read_from_file(fitxiz, LISTA_OBJETOS);
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
    case 9: /* <TAB> */
        siguiente_elemento_lista();
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
    mlist matriz_transformacion;

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
    camara_ptr = 0;
    obj_ptr = 0;
    camara_activa = 1;
    tipo_camara = CAMARA_PERSPECTIVA;
    ultimo_es_visible = 0;
    dibujar_no_visible = 1;
    dibujar_normales = 1;

    // TODO: Temporal, cargar una camara como un objeto. Buscar otra manera mas simple
    read_from_file("cam.obj", LISTA_CAMARAS);
    camara_ptr = (*sel_ptr);

    translacion(&matriz_transformacion, EJE_Z, DIR_ADELANTE, 300);
    aplicar_transformacion(&matriz_transformacion, SISTEMA_LOCAL);

    if (argc > 1)
        read_from_file(argv[1], LISTA_OBJETOS);
    else
    {
        read_from_file("r_falke.obj", LISTA_OBJETOS);
        translacion(&matriz_transformacion, EJE_X, DIR_ATRAS, 240);
        aplicar_transformacion(&matriz_transformacion, SISTEMA_LOCAL);

        read_from_file("x_wing.obj", LISTA_OBJETOS);
        translacion(&matriz_transformacion, EJE_Y, DIR_ATRAS, 20);
        aplicar_transformacion(&matriz_transformacion, SISTEMA_LOCAL);

        read_from_file("cam.obj", LISTA_OBJETOS);
        translacion(&matriz_transformacion, EJE_X, DIR_ADELANTE, 370);
        aplicar_transformacion(&matriz_transformacion, SISTEMA_LOCAL);

        read_from_file("box.obj", LISTA_OBJETOS);
        translacion(&matriz_transformacion, EJE_Z, DIR_ADELANTE, 200);
        translacion(&matriz_transformacion, EJE_X, DIR_ADELANTE, 200);
        aplicar_transformacion(&matriz_transformacion, SISTEMA_LOCAL);
    }

    mperspectiva_ptr = (mlist *)malloc(sizeof(mlist));
    calcular_mperspectiva(mperspectiva_ptr, CAMARA_CONFIG_NEAR, CAMARA_CONFIG_FAR, CAMARA_CONFIG_RIGHT, CAMARA_CONFIG_LEFT, CAMARA_CONFIG_TOP, CAMARA_CONFIG_BOTTOM);

    glutMainLoop();

    return 0;
}
