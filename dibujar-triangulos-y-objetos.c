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
#include "cargar-triangulo.h"

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
int denak;
int lineak;
int objektuak;
char aldaketa;
int ald_lokala;

char fitxiz[100];

void objektuari_aldaketa_sartu_ezk(double m[16])
{
}

void objektuari_aldaketa_sartu_esk(double m[16])
{
}

// TODO
// funtzio honek u eta v koordenatuei dagokien pointerra itzuli behar du.
// debe devolver el pointer correspondiente a las coordenadas u y v
unsigned char *color_textura(float u, float v)
{
    int desplazamendua;
    char *lag;

    int xp, yp;

    xp = u * dimx;
    yp = (1-v) * dimy; 

    desplazamendua = 1;
    lag = (unsigned char *)bufferra; // pixel on the left and top
    return (lag + dimx * yp * 3 + xp * 3); // Hay que multiplicar por 3 porque cada punto tiene (r,g,b)
    // return (lag + 3 * desplazamendua);
}

// TODO
// lerroa marrazten du, baina testuraren kodea egokitu behar da
// dibuja una linea pero hay que codificar la textura
void dibujar_linea_z(int linea, float c1x, float c1z, float c1u, float c1v, float c2x, float c2z, float c2u, float c2v)
{
    float xkoord, zkoord;
    float u, v;
    unsigned char r, g, b;
    unsigned char *colorv;

    glBegin(GL_POINTS);
    for (xkoord = c1x, zkoord = c1z, u = c1u, v = c1v; xkoord <= c2x; xkoord++)
    {
        // TODO
        // color_textura funtzioa ondo kodetu
        // programar de forma correcta la función color_textura
        colorv = color_textura(u, v);
        r = colorv[0];
        g = colorv[1];
        b = colorv[2];
        glColor3ub(r, g, b);
        glVertex3f(xkoord, linea, zkoord);
        // TODO
        // zkoord, u eta v berriak kalkulatu eskuineko puntuarentzat
        // calcular zkoord, u y v del siguiente pixel
    }
    glEnd();
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
    pptr->x = p.x;
    pptr->y = p.y;
    pptr->z = p.z;
    pptr->u = p.u;
    pptr->v = p.v;
}

void obtener_punto_corte(punto *pcorteptr, punto *psupptr, punto *pinfptr, int altura)
{

    /*
        /\y ---> /\y'
        /\x ---> /\x'

        /\y = pinf->y - psup->y
        /\y' = pcorte->y - psup->y

        /\x = pinf->x - psup->x
        /\x' = pcorte->x - psup->x

        Queremos calcula pcorte->x
        Hacemos la regla de tres para conseguir /\x'
        Y extraemos pcorte->x de la ecuación
    */

    pcorteptr->y = altura;

    int diffY = pinfptr->y - psupptr->y;    // deltaY
    int diffYp = pcorteptr->y - psupptr->y; // deltaY'

    int diffX = pinfptr->x - psupptr->x; // deltaX
    // int diffXp = pcorteptr->x - psupptr->x; // deltaX'
    int diffXp = (diffX * diffYp) / diffY;

    pcorteptr->x = diffXp + psupptr->x;

    /*
        Hemos calculado (x,y)
        Ahora vamos a calcular (z,u,v) usando el mismo metodo

        /\y ---> /\y'
        /\z ---> /\z'

        /\y ---> /\y'
        /\u ---> /\u'

        /\y ---> /\y'
        /\v ---> /\v'

    */

    int diffZ = pinfptr->z - psupptr->z;
    // int diffZp = pcorteptr->z - psupptr->z;
    int diffZp = (diffZ * diffYp) / diffY;
    pcorteptr->z = diffZp + psupptr->z;

    int diffU = pinfptr->u - psupptr->u;
    // int diffUp = pcorteptr->u - psupptr->u;
    int diffUp = (diffU * diffYp) / diffY;
    pcorteptr->u = diffU + psupptr->u;

    int diffV = pinfptr->v - psupptr->v;
    // int diffVp = pcorteptr->v - psupptr->v;
    int diffVp = (diffV * diffYp) / diffY;
    pcorteptr->v = diffV + psupptr->v;
}

int max(int *listptr, int length)
{
    int i;
    int max = listptr[0];

    for (i = 0; i < length; i++)
    {
        if (listptr[i] > max)
            max = *(listptr + i);
    }

    return max;
}

int min(int *listptr, int length)
{
    int i;
    int min = listptr[0];

    for (i = 0; i < length; i++)
    {
        if (listptr[i] < min)
            min = *(listptr + i);
    }

    return min;
}

void calcular_maximo_minimo(hiruki *tptr, punto *pmin, punto *pmax)
{
    int listax[3];
    int listay[3];
    int listaz[3];

    listax[0] = tptr->p1.x;
    listax[1] = tptr->p2.x;
    listax[2] = tptr->p3.x;

    listay[0] = tptr->p1.y;
    listay[1] = tptr->p2.y;
    listay[2] = tptr->p3.y;

    listaz[0] = tptr->p1.z;
    listaz[1] = tptr->p2.z;
    listaz[2] = tptr->p3.z;

    pmin->x = min(listax, 3);
    pmin->y = min(listay, 3);
    pmin->z = min(listaz, 3);

    pmax->x = max(listax, 3);
    pmax->y = max(listay, 3);
    pmax->z = max(listaz, 3);
}

void dibujar_triangulo(triobj *optr, int i)
{
    hiruki *tptr;

    punto *pgoiptr, *pbeheptr, *perdiptr;
    float x1, h1, z1, u1, v1, x2, h2, z2, u2, v2, x3, h3, z3, u3, v3;
    float c1x, c1z, c1u, c1v, c2x, c2z, c2u, c2v;
    int linea;
    float cambio1, cambio1z, cambio1u, cambio1v, cambio2, cambio2z, cambio2u, cambio2v;
    punto p1, p2, p3;

    float t, s, j;
    float cambiot, cambios, cambioj;

    // punto pcorte1, pcorte2;
    punto pmin, pmax;

    punto pcalculado, pcorte1, pcorte2;
    float alfa, beta, gamma;

    punto *pcortemayor, *pcortemenor;

    unsigned char r, g, b;
    unsigned char *colorv;

    if (i >= optr->num_triangles)
        return;
    tptr = optr->triptr + i;
    mxp(&p1, optr->mptr->m, tptr->p1);
    mxp(&p2, optr->mptr->m, tptr->p2);
    mxp(&p3, optr->mptr->m, tptr->p3);
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
    // hemen azpikoa kendu eta triangelua testurarekin marrazten duen kodea sartu.
    // lo que sigue aqui hay que sustituir por el código adecuado que dibuja el triangulo con textura

    // Calcular Psup, Pmed, Pinf

    if (tptr->p1.y > tptr->p2.y)
    {
        pgoiptr = &tptr->p1;  // Psup <- P1
        pbeheptr = &tptr->p2; // Pinf <- P2
    }
    else
    {
        pgoiptr = &tptr->p2;  // Psup <- P2
        pbeheptr = &tptr->p1; // Pinf <- P1
    }

    if (tptr->p3.y > pgoiptr->y)
    {
        perdiptr = pgoiptr;  // Pmed <- Psup
        pgoiptr = &tptr->p3; // Psup <- P3
    }
    else if (tptr->p3.y < pbeheptr->y)
    {
        perdiptr = pbeheptr;  // Pmed <- Pinf
        pbeheptr = &tptr->p3; // Psup <- P3
    }
    else
        perdiptr = &tptr->p3; // Pmed <- P3

    if (pgoiptr->y - perdiptr->y == 0)
        cambiot = 1;
    else
        cambiot = 1 / (pgoiptr->y - perdiptr->y);

    if (pgoiptr->y - pbeheptr->y == 0)
        cambios = 1;
    else
        cambios = 1 / (pgoiptr->y - pbeheptr->y);

    // De A -> B (t)  y  A -> C (s)
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

        if (pcorte1.x > pcorte2.x)
        {
            pcortemayor = &pcorte1;
            pcortemenor = &pcorte2;
        }
        else
        {
            pcortemayor = &pcorte2;
            pcortemenor = &pcorte1;
        }

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

    // Dibujar parte inferior
    // Seguimos de la s anterior
    // De B -> C (t)  y  A -> C (s)

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

        if (pcorte1.x > pcorte2.x)
        {
            pcortemayor = &pcorte1;
            pcortemenor = &pcorte2;
        }
        else
        {
            pcortemayor = &pcorte2;
            pcortemenor = &pcorte1;
        }

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
        optr->mptr->m[0] = 1.0;
        optr->mptr->m[5] = 1.0;
        optr->mptr->m[10] = 1.0;
        optr->mptr->m[15] = 1.0;
        optr->mptr->hptr = 0;
        // printf("objektu zerrendara doa informazioa...\n");
        optr->hptr = foptr;
        foptr = optr;
        sel_ptr = optr;
    }
    printf("datuak irakurrita\nLecura finalizada\n");
}

void x_aldaketa(int dir)
{
}

void y_aldaketa(int dir)
{
}

void z_aldaketa(int dir)
{
}

void undo()
{
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
    case 'g':
        if (ald_lokala == 1)
            ald_lokala = 0;
        else
            ald_lokala = 1;
        break;
    case 'x':
        x_aldaketa(1);
        break;
    case 'y':
        y_aldaketa(1);
        break;
    case 'z':
        z_aldaketa(1);
        break;
    case 'X':
        x_aldaketa(0);
        break;
    case 'Y':
        y_aldaketa(0);
        break;
    case 'Z':
        z_aldaketa(0);
        break;
    case 'u':
        undo();
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
    denak = 0;
    lineak = 0;
    objektuak = 0;
    foptr = 0;
    sel_ptr = 0;
    aldaketa = 'r';
    ald_lokala = 1;
    if (argc > 1)
        read_from_file(argv[1]);
    else
        read_from_file("adibideak.txt");
    glutMainLoop();

    return 0;
}
