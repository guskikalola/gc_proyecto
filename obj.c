#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
// #include "definitions.h"
#include "obj.h"

#define MAXLINE 200

void calcular_vnormal_cara(object3d *optr, int i)
{
    face *faceptr;
    faceptr = optr->face_table + i;

    // p1p2 /\ p1p3 = normal
    double p1p2[3] = {optr->vertex_table[faceptr->vertex_ind_table[1]].coord.x - optr->vertex_table[faceptr->vertex_ind_table[0]].coord.x, optr->vertex_table[faceptr->vertex_ind_table[1]].coord.y - optr->vertex_table[faceptr->vertex_ind_table[0]].coord.y, optr->vertex_table[faceptr->vertex_ind_table[1]].coord.z - optr->vertex_table[faceptr->vertex_ind_table[0]].coord.z}; // p2 - p1
    double p1p3[3] = {optr->vertex_table[faceptr->vertex_ind_table[2]].coord.x - optr->vertex_table[faceptr->vertex_ind_table[0]].coord.x, optr->vertex_table[faceptr->vertex_ind_table[2]].coord.y - optr->vertex_table[faceptr->vertex_ind_table[0]].coord.y, optr->vertex_table[faceptr->vertex_ind_table[2]].coord.z - optr->vertex_table[faceptr->vertex_ind_table[0]].coord.z}; // p3 - p1
    // double p1p3[3] = {faceptr->p3.x - faceptr->p1.x, faceptr->p3.y - faceptr->p1.y, faceptr->p3.z - faceptr->p1.z}; // p3 - p1
    double p1p2_p1p3[3] = {p1p2[1] * p1p3[2] - p1p2[2] * p1p3[1], -(p1p2[0] * p1p3[2] - p1p2[2] * p1p3[0]), p1p2[0] * p1p3[1] - p1p2[1] * p1p3[0]};

    // double mod_vp_zc = sqrt(pow(vp_zc[0], 2) + pow(vp_zc[1], 2) + pow(vp_zc[2], 2));
    double mod_p1p2_p1p3 = sqrt(pow(p1p2_p1p3[0], 2) + pow(p1p2_p1p3[1], 2) + pow(p1p2_p1p3[2], 2));

    faceptr->N[0] = p1p2_p1p3[0] / mod_p1p2_p1p3;
    faceptr->N[1] = p1p2_p1p3[1] / mod_p1p2_p1p3;
    faceptr->N[2] = p1p2_p1p3[2] / mod_p1p2_p1p3;
}

/*
 * Auxiliar function to process each line of the file
 */
static int sreadint(char *lerroa, int *zenbakiak)
{
    char *s = lerroa;
    int i, zbk, kont = 0;

    while (sscanf(s, " %d%n", &zbk, &i) > 0)
    {
        s += i;
        zenbakiak[kont++] = zbk;
    }
    return (kont);
}

static int sreadint2(char *lerroa, int *zenbakiak)
{
    char *s = lerroa;
    int i, zbk, kont = 0;

    while (sscanf(s, " %d%n", &zbk, &i) > 0)
    {
        s += i;
        while ((*s != ' ') && (*s != '\0'))
            s++; // jump vector normal information
        zenbakiak[kont++] = zbk;
    }
    printf("%d numbers in the line\n", kont);
    return (kont);
}
/**
 * @brief Function to read wavefront files (*.obj)
 * @param file_name Path of the file to be read
 * @param object_ptr Pointer of the object3d type structure where the data will be stored
 * @return Result of the reading: 0=Read ok, 1=File not found, 2=Invalid file, 3=Empty file
 */
int read_wavefront(char *file_name, object3d *object_ptr)
{
    vertex *vertex_table;
    face *face_table;
    int num_vertices = -1, num_faces = -1, count_vertices = 0, count_faces = 0;
    FILE *obj_file;
    char line[MAXLINE], line_1[MAXLINE], aux[45];
    int k;
    int i, j;
    int r, g, b;
    int values[MAXLINE];

    int tiene_color = 0;

    int num_vertices_face = 0;
    int indexx_face = 0;

    face *fptr;
    vertex *vptr;

    /*
     * The function reads twice the file. In the first read the number of
     * vertices and faces is obtained. Then, memory is allocated for each
     * of them and in the second read the actual information is read and
     * loaded. Finally, the object structure is created
     */
    if ((obj_file = fopen(file_name, "r")) == NULL)
        return (1);
    while (fscanf(obj_file, "\n%[^\n]", line) > 0)
    {
        i = 0;
        while (line[i] == ' ')
            i++;
        if ((line[0] == '#') && ((int)strlen(line) > 5))
        {
            i += 2;
            j = 0;
            // it is posible a line of the form "# number vertices" where "number" is a number
            // it is posible a line of the form "# number elements" where "number" is a number
            // it is posible a line of the form "# color r g b" where "r", "g" and "b" are numbers <256
            while (line[i] != ' ')
                line_1[j++] = line[i++];
            i++;
            line_1[j] = '\0';
            j = 0;
            if ((strcmp(line_1, "color") == 0) || (strcmp(line_1, "colour") == 0))
            {
                k = sscanf(line + i, "%d%d%d", &r, &g, &b);
                if (k == 3)
                {
                    object_ptr->rgb.r = r;
                    object_ptr->rgb.g = g;
                    object_ptr->rgb.b = b;

                    tiene_color = 1;
                }
            }
            else
            {
                while ((line[i] != ' ') && (line[i] != '\0'))
                    aux[j++] = line[i++];
                aux[j] = 0;
                if (strcmp(aux, "vertices") == 0)
                    num_vertices = atoi(line_1);
                if (strncmp(aux, "elements", 7) == 0)
                    num_faces = atoi(line_1);
            }
        }
        else
        {
            if (strlen(line) > 6)
            {
                if (line[i] == 'f' && line[i + 1] == ' ')
                    count_faces++;
                else if (line[i] == 'v' && line[i + 1] == ' ')
                    count_vertices++;
            }
        }
    }
    fclose(obj_file);
    printf("1 pasada: num vert = %d (%d), num faces = %d(%d) \n", num_vertices, count_vertices, num_faces, count_faces);
    if ((num_vertices != -1 && num_vertices != count_vertices) || (num_faces != -1 && num_faces != count_faces))
    {
        printf("WARNING: full file format: (%s)\n", file_name);
        // return (2);
    }
    if (num_vertices == 0 || count_vertices == 0)
    {
        printf("No vertex found: (%s)\n", file_name);
        return (3);
    }
    if (num_faces == 0 || count_faces == 0)
    {
        printf("No faces found: (%s)\n", file_name);
        return (3);
    }

    num_vertices = count_vertices;
    num_faces = count_faces;

    vertex_table = (vertex *)malloc(num_vertices * sizeof(vertex));
    face_table = (face *)malloc(num_faces * sizeof(face) * 2); // * 2 porque vamos a dividir los cuadrados en triangulos

    obj_file = fopen(file_name, "r");
    k = 0;
    j = 0;
    indexx_face = 0;

    for (i = 0; i < num_vertices; i++)
        vertex_table[i].num_faces = 0;

    while (fscanf(obj_file, "\n%[^\n]", line) > 0)
    {
        switch (line[0])
        {
        case 'v':
            if (line[1] == ' ') // vn not interested
            {
                sscanf(line + 2, "%lf%lf%lf", &(vertex_table[k].coord.x),
                       &(vertex_table[k].coord.y), &(vertex_table[k].coord.z));
                k++;
            }
            break;

        case 'f':
            if (line[1] == ' ') // fn not interested
            {
                for (i = 2; i <= (int)strlen(line); i++)
                    line_1[i - 2] = line[i];
                line_1[i - 2] = '\0';
                num_vertices_face = sreadint2(line_1, values);

                // Dividir poligono en triangulos
                // https://stackoverflow.com/a/23724231
                // 0 (i) (i + 1)  [for i in 1..(n - 2)]
                for (i = 1; i <= num_vertices_face-2; i++)
                {
                    face_table[indexx_face].num_vertices = 3; // Vamos a divirlo en triangulos, 3 vertices

                    face_table[indexx_face].vertex_ind_table = (int *)malloc(face_table[indexx_face].num_vertices * sizeof(int));
                    face_table[indexx_face].vertex_ind_table[0] = values[0] - 1;
                    face_table[indexx_face].vertex_ind_table[1] = values[i] - 1;
                    face_table[indexx_face].vertex_ind_table[2] = values[i+1] - 1;

                    // printf("f %d vertices\n",face_table[j].num_vertices);
                    for (j = 0; j < face_table[indexx_face].num_vertices; j++)
                    {
                        // face_table[indexx_face].vertex_ind_table[j] = values[j] - 1;
                        // printf(" %d ",values[i] - 1);
                        vertex_table[face_table[indexx_face].vertex_ind_table[j]].num_faces++;
                    }

                    indexx_face++;
                }

                // printf("\n");
                // }

                j++;
            }
            break;
        }
    }

    num_faces = indexx_face;

    fclose(obj_file);

    printf("2 pasada\n");

    /*
     * Information read is introduced in the structure */
    object_ptr->vertex_table = vertex_table;
    object_ptr->face_table = face_table;
    object_ptr->num_vertices = num_vertices;
    object_ptr->num_faces = num_faces;

    /*
     * The maximum and minimum coordinates are obtained **/
    object_ptr->max.x = object_ptr->vertex_table[0].coord.x;
    object_ptr->max.y = object_ptr->vertex_table[0].coord.y;
    object_ptr->max.z = object_ptr->vertex_table[0].coord.z;
    object_ptr->min.x = object_ptr->vertex_table[0].coord.x;
    object_ptr->min.y = object_ptr->vertex_table[0].coord.y;
    object_ptr->min.z = object_ptr->vertex_table[0].coord.z;

    for (i = 1; i < object_ptr->num_vertices; i++)
    {
        if (object_ptr->vertex_table[i].coord.x < object_ptr->min.x)
            object_ptr->min.x = object_ptr->vertex_table[i].coord.x;

        if (object_ptr->vertex_table[i].coord.y < object_ptr->min.y)
            object_ptr->min.y = object_ptr->vertex_table[i].coord.y;

        if (object_ptr->vertex_table[i].coord.z < object_ptr->min.z)
            object_ptr->min.z = object_ptr->vertex_table[i].coord.z;

        if (object_ptr->vertex_table[i].coord.x > object_ptr->max.x)
            object_ptr->max.x = object_ptr->vertex_table[i].coord.x;

        if (object_ptr->vertex_table[i].coord.y > object_ptr->max.y)
            object_ptr->max.y = object_ptr->vertex_table[i].coord.y;

        if (object_ptr->vertex_table[i].coord.z > object_ptr->max.z)
            object_ptr->max.z = object_ptr->vertex_table[i].coord.z;
    }

    printf("Calculando normales...\n");

    for (i = 0; i < num_faces; i++)
    {

        // Calcular vector normal de la cara
        calcular_vnormal_cara(object_ptr, i);

        fptr = object_ptr->face_table + i;

        // Sumar normal de la cara a los vertices que usa
        for (j = 0; j < fptr->num_vertices; j++)
        {
            vptr = object_ptr->vertex_table + fptr->vertex_ind_table[j];
            // printf("----(S)-----\n");
            // printf("fptr->N[0] = %f, fptr->N[1] = %f, fptr->N[2] = %f\n", fptr->N[0], fptr->N[1], fptr->N[2]);
            // printf("vptr->N[0] = %f, vptr->N[1] = %f, vptr->N[2] = %f\n", vptr->N[0], vptr->N[1], vptr->N[2]);
            vptr->N[0] = vptr->N[0] + fptr->N[0];
            vptr->N[1] = vptr->N[1] + fptr->N[1];
            vptr->N[2] = vptr->N[2] + fptr->N[2];
            // printf("fptr->N[0] = %f, fptr->N[1] = %f, fptr->N[2] = %f\n", fptr->N[0], fptr->N[1], fptr->N[2]);
            // printf("vptr->N[0] = %f, vptr->N[1] = %f, vptr->N[2] = %f\n", vptr->N[0], vptr->N[1], vptr->N[2]);
            // printf("----(E)-----\n");
        }
    }

    // Normalizar vectores "normales" ( aun no lo son ) de los vertices
    for (i = 0; i < num_vertices; i++)
    {
        vptr = object_ptr->vertex_table + i;

        double mod_n = sqrt(pow(vptr->N[0], 2) + pow(vptr->N[1], 2) + pow(vptr->N[2], 2));

        vptr->N[0] = vptr->N[0] / mod_n;
        vptr->N[1] = vptr->N[1] / mod_n;
        vptr->N[2] = vptr->N[2] / mod_n;

        // printf("mod_n=%f\n", mod_n);
    }

    if (tiene_color == 0)
    {
        object_ptr->rgb.r = 0;
        object_ptr->rgb.g = 0;
        object_ptr->rgb.b = 0;
    }

    return (0);
}
